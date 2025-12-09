# Publishing Guide for pioud_soap

This document provides instructions on how to build and publish the `pioud_soap` package to PyPI.

## Prerequisites

Make sure you have the following tools installed:

```bash
pip install build twine
```

## Building the Package

To build the package, run the following command from the root directory of the project:

```bash
python -m build
```

This will create both source distribution (`.tar.gz`) and wheel (`.whl`) files in the `dist/` directory.

## Testing the Package Locally

Before publishing, you can test the package locally by installing it in development mode:

```bash
pip install -e .
```

Or by installing the built wheel:

```bash
pip install dist/pioud_soap-0.1.0-py3-none-any.whl
```

## Publishing to TestPyPI (Recommended for Testing)

Before publishing to the main PyPI repository, it's a good practice to test on TestPyPI:

```bash
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

You can then install from TestPyPI to verify everything works:

```bash
pip install --index-url https://test.pypi.org/simple/ pioud_soap
```

## Publishing to PyPI

Once you're confident that the package works correctly, you can publish it to the main PyPI repository:

```bash
twine upload dist/*
```

You'll need to provide your PyPI username and password.

## Updating the Package

To update the package:

1. Update the version number in `pioud_soap/__init__.py` and `setup.py`
2. Rebuild the package: `python -m build`
3. Upload the new version: `twine upload dist/*`

## Creating a GitHub Release

After publishing to PyPI, it's a good practice to create a release on GitHub:

1. Tag the commit: `git tag v0.1.0`
2. Push the tag: `git push origin v0.1.0`
3. Create a release on GitHub with release notes