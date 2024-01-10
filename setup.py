from setuptools import setup, find_packages

setup(
      name='basefly',
      version='1.0.0',
      url='https://github.com/gudeqing/basefly',
      license='MIT',
      author='gudeqing',
      author_email='822466659@qq.com',
      description='Easy to create workflow for heavy computing',
      packages=find_packages(include=['basefly']),
      long_description=open('README.md', encoding='utf-8').read(),
      zip_safe=False,
      classifiers=[
            "Development Status :: 1 - Alpha",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: bio informatics"
        ],
      install_requires=["argparse>=1.1"],
      setup_requires=[],
      python_requires=">=3.7"
)

# python setup.py sdist bdist_wheel
