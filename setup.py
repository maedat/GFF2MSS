from setuptools import setup, find_packages

setup(
    name="gff2mss",
    version="4.1.1.1",
    author="Taro Maeda",
    author_email="taromaedaj@gmail.com",
    description="Convert GFF3 gene models to MSS format for DDBJ submission",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/maedat/GFF2MSS",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "biopython",
        "BCBio",
        "gffpandas"
    ],
    entry_points={
        "console_scripts": [
            "gff2mss = GFF2MSS_Improved:main"  # スクリプトをコマンドラインツールとして登録
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)