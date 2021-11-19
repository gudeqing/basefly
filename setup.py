from setuptools import setup, find_packages

setup(
      name='basefly',
      version='1.0.0',
      url='https://github.com/gudeqing/nestcmd',
      license='MIT',
      author='gudeqing',
      author_email='822466659@qq.com',
      description='building workflow for heavy computing',
      packages=find_packages(exclude=['tests', 'docs']),
      long_description=open('README.md').read(),
      zip_safe=False,
      classifiers=[
            "Development Status :: 1 - Alpha",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: bio informatics"
        ],
      install_requires=["argparse>=1.1.0", "munch>=2.5.0"],
      setup_requires=[],
)
