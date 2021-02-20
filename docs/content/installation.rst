############
Installation
############

====================
PyPi Install
====================
| 1)  Install xpressplot and associated dependencies via pip:

.. code-block:: shell

  $ pip install xpressplot

====================
Conda Install
====================
| This feature is not yet available...
|
| 1)  Install xpressplot and associated dependencies via conda:

.. code-block:: shell

  $ conda install -y -c bioconda xpressplot

======================
Manual install
======================

| 1)  Or download xpressplot manually:

.. code-block:: shell

  $ git clone https://github.com/XPRESSyourself/xpressplot.git
  $ cd xpressplot
  $ python setup.py install

| 2)  Or, to download specific version:

.. code-block:: shell

  $ tag='v0.0.1-beta'
  $ wget https://github.com/XPRESSyourself/xpressplot/archive/$tag.zip
  $ unzip xpressplot-${tag:1}.zip
  $ mv xpressplot-${tag:1} xpressplot
  $ cd xpressplot
  $ python setup.py install

| 3)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing xpressplot script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/xpressplot-0.0.1b0-py3.6.egg
  Processing dependencies for xpressplot==0.0.1b0
  Finished processing dependencies for xpressplot==0.0.1b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile
