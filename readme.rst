
gtex
====

A small collection of scripts to process and format GTEx__ eQTL datasets.

.. __: https://gtexportal.org/home/index.html

The file ``src/process-brain-eqtls.sh`` contains a sample processing pipeline for 
brain tissue eQTLs.


Usage
-----

Retrieve GTEx eQTLs v.7, the GTEx lookup table, and annotations.
eQTLs and metadata are decompressed and stored in ``data/``. 

.. code:: bash

    $ ./src/get-gtex-eqtls.sh

Process GTEx eQTLs into a tab delimited format with tissue annotations and dbSNP
reference identifiers:

.. code:: bash

    $ python src/process_gtex_eqtls.py -d -a data/gtex-v7-annotations.tsv data/gtex-eqtls data/gtex-lookup-table.tsv.gz processed-gtex-eqtls.tsv

GTEx uses dbSNP v.147 to annotate eQTL variants.
If you need a more recent version of dbSNP (e.g., v.150) you can use the merge
scripts to update SNP identifiers.

First, retrieve the merge table from NCBI:

.. code:: bash

    $ ./src/get-merge-table.sh --v150

Then use the merge script to update GTEx SNP IDs:

.. code:: bash

    $ ./src/merge_snp_ids.py processed-gtex-eqtls.tsv data/dbsnp-merge-table.bcp.gz processed-gtex-eqtls-v150.tsv


Requirements and installation
-----------------------------

The following software and packages are required:

- Python 2.7/3.6/3.7
- NumPy__
- Pandas__

.. __: https://www.numpy.org/
.. __: https://pandas.pydata.org/

To install, use pip:

.. code:: bash

    $ pip install -r requirements.txt

