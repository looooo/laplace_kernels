from distutils.core import setup, Extension


ext2D = Extension("laplace/laplaceKern2D",

  )
ext3D = Extension()
extVtkWriter = Extension()

setup(name='ppm._ppm',
      version='0.0.1',
      author='Loooo',
      requires='eigen',
      author_email='sppedflyer@gmail.com',
      url="https://github.com/looooo/panelmethod",
      description='Wrap PM using pybind11',
      packages=files,
      ext_modules=[Extension('ppm._ppm',
                   sources=src,
                   include_dirs=include_dirs,
                   extra_compile_args=['-std=c++11', '-fopenmp'],
                   extra_link_args=extra_link_args)])
