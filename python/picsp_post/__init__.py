"""picsp_post: post-processing for PICSP3D output.

Quick start::

    from picsp_post import Run
    run = Run("output_dir")
    run.history()                    # numpy structured array
    grid = run.fields(-1)            # pyvista ImageData of the last dump
    phi = run.field_array(-1, "phi") # (nx+1, ny+1, nz+1) numpy array

Command line: ``python -m picsp_post --help``
"""
from .io import Run  # noqa: F401
from . import style  # noqa: F401

__version__ = "3.0.0"
