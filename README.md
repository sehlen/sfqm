Finite Quadratic Modules and Simple Lattices
============================================

Implementation of several algorithms, tools and supplementary material
for the article "Lattices with many Borcherds products" 
by Jan H. Bruinier, S. Ehlen and E. Freitag.

Installation
============

-  A) On your own computer which has sage and git installed:
   We expect that a current version of sage is installed on your system (see dependencies).

   Go to the directory where you would like to download the repository to.

   You should be able to simply run:
   ```
    $ git clone https://github.com/sehlen/sfqm.git
    $ sage setup.py install
   ```

   Then start sage and try the following:
   ```
   sage: from sfqm.fqm.genus_symbol import *
   sage: s = GenusSymbol('3^+1')
   sage: s
   Genus symbol 3^+1
   sage: s.dimension_cusp_forms(3)
   0
   ```

-  B) In the sage math cloud (SMC), http://cloud.sagemath.com:
   This is a nice possibility to try out the programs
   if you don't have sage installed.

   You only have to modify A) a little bit.
   Instead of the first command you can (in a project of your choice)
   click on the "+ New" button whcih you use to create a new file.

   Then you insert the url ```https://github.com/sehlen/sfqm.git```
   in the input field (filename) and click on "From Internet".

   Now open a new terminal ("+New" and then "Terminal") and type:

   ```
    $ cd sfqm
    $ sage setup.py install --user
   ```

   Now you can test the installation as in A),
   either on the command line or in a sage worksheet 
   (you may have to restart the worksheet server in your project).
   

Dependencies
============

- sage: A current version of sage 
(we tested the implementation with sage 6.1.1 and sage 6.2)

- git: if you do not have git installed, you can download the tarball version.

The library dependencies have been included into the release files.

The following modules are required and can also be obtained
separately from my (sehlen) psage fork (use the devel branch).

- psage.modules.finite_quadratic_module
- psage.modform.weilrep_tools
- psage.modules.weil_invariants
