############
Installation
############

====================
PyPi Install
====================
| 1)  Install XPRESStools and associated dependencies via pip:

.. code-block:: shell

  $ pip install xpresstools

====================
Conda Install
====================
| This feature is not yet available...
|
| 1)  Install XPRESStools and associated dependencies via conda:

.. code-block:: shell

  $ conda install -y -c bioconda xpresstools

======================
Manual install
======================

| 1)  Or download XPRESStools manually:

.. code-block:: shell

  $ git clone https://github.com/XPRESSyourself/xpresstools.git
  $ cd XPRESStools
  $ cd python setup.py install

| 2)  Or, to download specific version:

.. code-block:: shell

  $ tag='v0.0.1-beta'
  $ wget https://github.com/XPRESSyourself/XPRESStools/archive/$tag.zip
  $ unzip xpresstools-${tag:1}.zip
  $ mv xpresstools-${tag:1} xpresstools
  $ cd xpresstools
  $ cd python setup.py install

| 3)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing xpresstools script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/xpresstools-0.0.1b0-py3.6.egg
  Processing dependencies for XPRESStools==0.0.1b0
  Finished processing dependencies for XPRESStools==0.0.1b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile
