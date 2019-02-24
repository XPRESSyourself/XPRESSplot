############
Installation
############

1)  You may need to manually point your system to python3. You can check this by typing :data:`python` in the command line and seeing if it is running python3. At this time, RiboPipe only works when Python3 is set as the default. If it is not, do the following:

.. code-block:: shell

  $ echo "alias python="python3" >> ~/.bash_profile

2)  Download `Conda <https://www.anaconda.com/download/>`_, a package manager, for your operating system. Double click the `.pkg` file if on MacOS, the `.exe` file on Windows, or follow these `instructions <https://conda.io/docs/user-guide/install/linux.html#install-linux-silent>`_ on Linux.

3)  Install XPRESStools and associated dependencies via conda:

.. code-block:: shell

  $ conda install -y -c bioconda xpresstools

4)  Or download XPRESStools manually:

.. code-block:: shell

  $ git clone https://github.com/XPRESSyourself/xpresstools.git
  $ cd XPRESStools
  $ cd python setup.py install

4)  Or, to download specific version:

.. code-block:: shell

  $ tag='v0.0.1-beta'
  $ wget https://github.com/XPRESSyourself/XPRESStools/archive/$tag.zip
  $ unzip xpresstools-${tag:1}.zip
  $ mv xpresstools-${tag:1} xpresstools
  $ cd xpresstools
  $ cd python setup.py install

9)  At the end of the installation instructions, an installation location will be given. Add this to your $PATH:

.. code-block:: shell

  ...
  Installing xpresstools script to /Users/$USERNAME/anaconda3/bin

  Installed /Users/$USERNAME/anaconda3/lib/python3.6/site-packages/xpresstools-0.0.1b0-py3.6.egg
  Processing dependencies for XPRESStools==0.0.1b0
  Finished processing dependencies for XPRESStools==0.0.1b0

  $ echo "export PATH='/Users/$USERNAME/anaconda3/bin:$PATH' >> ~/.bash_profile
