from setuptools import setup, find_packages

setup(
    name="transmorph",
    version="0.0.1", 
    packages=find_packages(),  
    install_requires=[  

        'pandas'
        'pysam'
    ],
    author="Yujia Feng",
    author_email="yfeng80@jh.edu",
    description="A tool that converts a BAM file with read alignments in genomic coordinates into transcriptomic coordinates using a given transcript annotation file",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yujiafeng8888/Transmorph",  
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.12',  
)
