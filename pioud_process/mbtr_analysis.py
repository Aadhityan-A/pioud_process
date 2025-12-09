import numpy as np
from ase.io import read, write
from ase.data import atomic_numbers  # Import atomic_numbers dictionary from ASE
from scipy.spatial.distance import pdist, squareform
from dscribe.descriptors import MBTR
from skmatter.sample_selection import FPS
from scipy.linalg import svd
from numpy.linalg import pinv
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def mbtr_analysis_from_pioudop(path, parsed_MBTR_dict):
    """
    Analyze structures using MBTR descriptors and select distinct structures using FPS and CUR.
    
    Parameters:
    -----------
    path : str
        Path to the directory containing input files and where output files will be written.
    parsed_MBTR_dict : dict
        Dictionary containing MBTR parameters and analysis settings.
    
    Returns:
    --------
    dict or None
        Dictionary containing output file paths and selected indices.
    """
    # Extract parameters from input dictionary
    n_to_select = int(parsed_MBTR_dict.get("n_to_select"))
    input_file_name = parsed_MBTR_dict.get("input_file_name")
    output_file_name = parsed_MBTR_dict.get("output_file_name")
    no_of_jobs = parsed_MBTR_dict.get("no_of_jobs")
    periodic = parsed_MBTR_dict.get("periodic", True)
    
    # Determine selection method: 'FPS', 'CUR', or 'CURFPS'
    method = parsed_MBTR_dict.get("method", "CURFPS").upper()
    if method not in ("FPS", "CUR", "CURFPS"):
        raise ValueError(
            f"Invalid method '{method}'. Choose 'FPS', 'CUR' or 'CURFPS'."
        )
    print(f"Selection method: {method}")

    # MBTR specific parameters with defaults
    # default_k2 = {
    #     "geometry": {"function": "inverse_distance"},
    #     "grid": {"min": 0.1, "max": 2, "n": 50, "sigma": 0.1}
    # }
    
    # default_k3 = {
    #     "geometry": {"function": "angle"},
    #     "grid": {"min": 0, "max": 180, "n": 50, "sigma": 5}
    # }
    
    # k2 = parsed_MBTR_dict.get("k2", default_k2)
    # k3 = parsed_MBTR_dict.get("k3", default_k3)
    
    # CUR parameters
    n_components = parsed_MBTR_dict.get("n_components", min(10, n_to_select))
    
    # Print important parameters
    print(f"Input parameters:")
    print(f"  Number of jobs: {no_of_jobs}")
    print(f"  Periodic: {periodic}")

    print(f"  n_components for CUR: {n_components}")
    print(f"path: {path}")

    # Reading the structure from file
    print(f"Reading structure from {input_file_name}...")
    structures = read(path+input_file_name, index=':')
    print(f"Number of structures read: {len(structures)}")

    # Collecting chemical species
    print("Collecting chemical species from structures...")
    species = set()
    for structure in structures:
        species.update(structure.get_chemical_symbols())
    species = list(species)
    print(f"Species found: {species}")

    # Extract MBTR parameters from parsed_MBTR_dict or use defaults
    geometry = {
        "function": parsed_MBTR_dict.get("geometry_function", "inverse_distance")
    }
    grid = {
        "min": parsed_MBTR_dict.get("grid_min", 0),
        "max": parsed_MBTR_dict.get("grid_max", 1),
        "n": parsed_MBTR_dict.get("grid_n", 100),
        "sigma": parsed_MBTR_dict.get("grid_sigma", 0.1)
    }
    weighting = {
        "function": parsed_MBTR_dict.get("weighting_function", "exp"),
        "scale": parsed_MBTR_dict.get("weighting_scale", 0.5),
        "threshold": parsed_MBTR_dict.get("weighting_threshold", 1e-3)
    }
    normalization = parsed_MBTR_dict.get("normalization", "none")

    # Setting up MBTR descriptor
    print("Setting up MBTR descriptor...")
    mbtr = MBTR(
        species=species,
        periodic=periodic,
        geometry=geometry,
        grid=grid,
        weighting=weighting,
        normalization=normalization,
    )
    # Print number of features
    num_features = mbtr.get_number_of_features()
    print(f"MBTR descriptor set up complete.")
    print(f"Number of features: {num_features}")

    # Generating MBTR feature vectors
    print("Generating MBTR feature vectors for the structures...")
    feature_vectors = mbtr.create(structures, n_jobs=no_of_jobs)
    print("MBTR feature vectors generated.")
    print(f"Feature shape: {feature_vectors.shape}")

    # ----
    # CUR decomposition (if requested)
    # ----
    if method in ("CUR", "CURFPS"):
        print("Performing complete CUR decomposition on feature matrix...")
        # 1) SVD
        U_svd, s_svd, Vt_svd = svd(feature_vectors, full_matrices=False)
        print(f"  SVD shapes — U: {U_svd.shape}, Σ: {s_svd.shape}, Vt: {Vt_svd.shape}")

        # 2) Column leverage (c = n_components)
        col_lev = np.sum(Vt_svd[:n_components, :]**2, axis=0)
        selected_columns = np.argsort(col_lev)[-n_components:]
        print(f"  Selected {len(selected_columns)} columns (C indices): {selected_columns}")

        # 3) Row leverage (r = n_to_select)
        row_lev = np.sum(U_svd[:, :n_components]**2, axis=1)
        selected_rows = np.argsort(row_lev)[-n_to_select:]
        print(f"  Selected {len(selected_rows)} rows (R indices): {selected_rows}")

        # 4) Build C, R, and W
        C = feature_vectors[:, selected_columns]
        R = feature_vectors[selected_rows, :]
        W = feature_vectors[np.ix_(selected_rows, selected_columns)]
        print(f"  W shape: {W.shape}")

        # 5) CUR U = pinv(W)
        U_cur = pinv(W)
        print(f"  CUR U matrix shape: {U_cur.shape}")

        # (opt) approximate back
        A_cur_approx = C.dot(U_cur).dot(R)
        print(f"  CUR approx shape: {A_cur_approx.shape}")

        # Write CUR‐selected structures
        cur_indices = selected_rows
        cur_structures = [structures[i] for i in cur_indices]
        cur_output_file = path + output_file_name
        write(cur_output_file, cur_structures)
        print(f"CUR‐selected structures written to {cur_output_file}")
        if len(cur_indices) != len(set(cur_indices)):
            print("ERROR: Duplicates in CUR indices.")
        else:
            print("No duplicates in CUR indices.")
    # --- end CUR ---

    # FPS selection method
    # ----
    # FPS selection (if requested)
    # ----
    if method in ("FPS", "CURFPS"):
        print("Running FPS selection for the structures...")

        # Choose features based on method
        if method == "CURFPS":
            print("Using CUR-reduced features for FPS...")
            features_for_fps = A_cur_approx
        else:
            print("Using direct MBTR feature vectors for FPS...")
            features_for_fps = feature_vectors

        fps_selector = FPS(
            initialize=0,
            n_to_select=n_to_select,
            progress_bar=True,
            full=False
        )

        # Running the FPS selection
        fps_selector.fit(features_for_fps)
        fps_indices = fps_selector.selected_idx_
        print(f"FPS selected indices: {fps_indices}")

        # Write FPS-selected structures to output file
        fps_structures = [structures[i] for i in fps_indices]
        fps_output_file = path + output_file_name
        write(fps_output_file, fps_structures)
        print(f"FPS-selected structures written to {fps_output_file}.")

        # Check for duplicates
        if len(fps_indices) != len(set(fps_indices)):
            print("ERROR: Duplicates found in FPS indices.")
        else:
            print("No duplicates in FPS indices.")
    # --- end FPS ---

    if method == "CURFPS":
        # Compare FPS and CUR selections
        overlap = set(fps_indices).intersection(set(cur_indices))
        print(f"Overlap between FPS and CUR selections: {len(overlap)} structures")
        print(f"Overlap percentage: {len(overlap)/n_to_select*100:.2f}%")

    # Build return dict
    results = {}
    if method in ("FPS", "CURFPS"):
        results["fps_output"] = fps_output_file
        results["fps_indices"] = fps_indices.tolist()
    if method in ("CUR", "CURFPS"):
        results["cur_output"] = cur_output_file
        results["cur_indices"] = cur_indices.tolist()
    if method == "CURFPS":
        overlap = set(fps_indices).intersection(cur_indices)
        results["overlap"] = {
            "indices": list(overlap),
            "percentage": len(overlap) / n_to_select * 100
        }
    # Plot MBTR 2D‐PCA and highlight FPS, CUR, and CURFPS selections
    try:
        print("Generating 2D PCA visualization of MBTR feature space...")
        
        # 1) Apply PCA to reduce MBTR features to 2D for visualization
        # PCA finds the principal components that capture maximum variance
        pca = PCA(n_components=2)
        mbtr_2d = pca.fit_transform(feature_vectors)
        
        # Print PCA statistics for scientific validation
        explained_variance_ratio = pca.explained_variance_ratio_
        cumulative_variance = np.sum(explained_variance_ratio)
        print(f"  PCA explained variance ratio: {explained_variance_ratio}")
        print(f"  Cumulative variance explained by 2 components: {cumulative_variance:.4f} ({cumulative_variance*100:.2f}%)")
        
        # Validate PCA transformation
        print(f"  Original feature space: {feature_vectors.shape[1]}D")
        print(f"  Reduced feature space: {mbtr_2d.shape[1]}D")
        
        # 2) Direct FPS on raw MBTR features (not PCA-reduced)
        fps_direct = FPS(initialize=0, n_to_select=n_to_select, full=False)
        fps_direct.fit(feature_vectors)  # Use original high-dimensional features
        fps_idx_direct = fps_direct.selected_idx_

        # Create plots based on selected method
        if method == "FPS":
            plt.figure(figsize=(10, 8))
            plt.scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                      c="lightgray", s=30, alpha=0.6, label="All structures")
            plt.scatter(mbtr_2d[fps_idx_direct, 0], mbtr_2d[fps_idx_direct, 1],
                      c="red", marker="o", s=100, edgecolors='black', linewidth=1,
                      label=f"FPS selected ({len(fps_idx_direct)})")
            
            plt.xlabel(f"PC1 ({explained_variance_ratio[0]:.3f} variance)")
            plt.ylabel(f"PC2 ({explained_variance_ratio[1]:.3f} variance)")
            plt.title(f"2D PCA of MBTR Features - FPS Selection\n"
                     f"Total variance explained: {cumulative_variance:.3f}")
            plt.legend(loc="best")
            plt.grid(True, alpha=0.3)
            
            plot_file = path + "fps_selection_plot.png"
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            plt.close()
            
        elif method == "CUR":
            plt.figure(figsize=(10, 8))
            plt.scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                      c="lightgray", s=30, alpha=0.6, label="All structures")
            plt.scatter(mbtr_2d[cur_indices, 0], mbtr_2d[cur_indices, 1],
                      c="blue", marker="s", s=100, edgecolors='black', linewidth=1,
                      label=f"CUR selected ({len(cur_indices)})")
            
            plt.xlabel(f"PC1 ({explained_variance_ratio[0]:.3f} variance)")
            plt.ylabel(f"PC2 ({explained_variance_ratio[1]:.3f} variance)")
            plt.title(f"2D PCA of MBTR Features - CUR Selection\n"
                     f"Total variance explained: {cumulative_variance:.3f}")
            plt.legend(loc="best")
            plt.grid(True, alpha=0.3)
            
            plot_file = path + "cur_selection_plot.png"
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            plt.close()
            
        else:  # CURFPS
            # FPS on CUR-approximated features
            fps_curfps = FPS(initialize=0, n_to_select=n_to_select, full=False)
            fps_curfps.fit(A_cur_approx)  # Use CUR-approximated features
            fps_idx_curfps = fps_curfps.selected_idx_
            
            # Individual plots
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            
            # FPS only
            axes[0,0].scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                           c="lightgray", s=30, alpha=0.6, label="All structures")
            axes[0,0].scatter(mbtr_2d[fps_idx_direct, 0], mbtr_2d[fps_idx_direct, 1],
                           c="red", marker="o", s=100, edgecolors='black', linewidth=1,
                           label=f"FPS ({len(fps_idx_direct)})")
            axes[0,0].set_title("FPS Selection")
            axes[0,0].legend()
            axes[0,0].grid(True, alpha=0.3)
            
            # CUR only
            axes[0,1].scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                           c="lightgray", s=30, alpha=0.6, label="All structures")
            axes[0,1].scatter(mbtr_2d[cur_indices, 0], mbtr_2d[cur_indices, 1],
                           c="blue", marker="s", s=100, edgecolors='black', linewidth=1,
                           label=f"CUR ({len(cur_indices)})")
            axes[0,1].set_title("CUR Selection")
            axes[0,1].legend()
            axes[0,1].grid(True, alpha=0.3)
            
            # CURFPS only
            axes[1,0].scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                           c="lightgray", s=30, alpha=0.6, label="All structures")
            axes[1,0].scatter(mbtr_2d[fps_idx_curfps, 0], mbtr_2d[fps_idx_curfps, 1],
                           c="green", marker="^", s=100, edgecolors='black', linewidth=1,
                           label=f"CURFPS ({len(fps_idx_curfps)})")
            axes[1,0].set_title("CURFPS Selection")
            axes[1,0].legend()
            axes[1,0].grid(True, alpha=0.3)
            
            # Combined comparison
            axes[1,1].scatter(mbtr_2d[:, 0], mbtr_2d[:, 1],
                           c="lightgray", s=20, alpha=0.4, label="All structures")
            axes[1,1].scatter(mbtr_2d[fps_idx_direct, 0], mbtr_2d[fps_idx_direct, 1],
                           c="red", marker="o", s=80, alpha=0.8, label="FPS")
            axes[1,1].scatter(mbtr_2d[cur_indices, 0], mbtr_2d[cur_indices, 1],
                           c="blue", marker="s", s=80, alpha=0.8, label="CUR")
            axes[1,1].scatter(mbtr_2d[fps_idx_curfps, 0], mbtr_2d[fps_idx_curfps, 1],
                           c="green", marker="^", s=80, alpha=0.8, label="CURFPS")
            axes[1,1].set_title("Combined Comparison")
            axes[1,1].legend()
            axes[1,1].grid(True, alpha=0.3)
            
            # Set common labels
            for ax in axes.flat:
                ax.set_xlabel(f"PC1 ({explained_variance_ratio[0]:.3f} variance)")
                ax.set_ylabel(f"PC2 ({explained_variance_ratio[1]:.3f} variance)")
            
            plt.suptitle(f"2D PCA of MBTR Features - Method Comparison\n"
                       f"Total variance explained: {cumulative_variance:.3f} "
                       f"({cumulative_variance*100:.1f}%)", fontsize=14)
            plt.tight_layout()
            
            plot_file = path + "curfps_combined_selection_plot.png"
            plt.savefig(plot_file, dpi=300, bbox_inches="tight")
            plt.close()

        print(f"PCA visualization saved to {plot_file}")
        
        # Add PCA information to results
        results["pca_info"] = {
            "explained_variance_ratio": explained_variance_ratio.tolist(),
            "cumulative_variance": float(cumulative_variance),
            "original_dimensions": int(feature_vectors.shape[1]),
            "reduced_dimensions": 2
        }
        results["selection_plot"] = plot_file
        
    except Exception as plot_e:
        print(f"Warning: failed to generate PCA visualization: {plot_e}")
        import traceback
        print(traceback.format_exc())

    return output_file_name

if __name__ == "__main__":
    path = "/home/aadhityan/Documents/GitHub/pioud_process/test_files_1/"
    parsed_MBTR_dict = {
        "periodic": True,
        "n_to_select": 2,
        "no_of_jobs": 3,
        "n_components": 2,
        "method": "CURFPS",
        "input_file_name": "ase_config_for_mace.xyz",
        "output_file_name": "ase_config_for_mace_1.xyz",
        # MBTR parameter keys
        "geometry_function": "inverse_distance",
        "grid_min": 0,
        "grid_max": 1,
        "grid_n": 100,
        "grid_sigma": 0.1,
        "weighting_function": "exp",
        "weighting_scale": 0.5,
        "weighting_threshold": 1e-3,
        "normalization": "none"
    }
    mbtr_analysis_from_pioudop(path, parsed_MBTR_dict)
