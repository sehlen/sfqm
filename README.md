Finite Quadratic Modules and Simple Lattices
============================================

Implementation of several algorithms and tools for the article 
"Lattices with many Borcherds products" 
by Jan H. Bruinier, S. Ehlen and E. Freitag.

Dependencies
============

- sage: A current version of sage 
(we tested the implementation with sage 6.1.1 and sage 6.2)

The library dependencies have been included as a git subtree.
So if you do not have psage installed or if you use the
main psage repository (wstein/psage), you can install the
required packages directly from this repository.

The following modules are required and can also be obtained
separately from my (sehlen) psage fork (use the devel branch).

- psage.modules.finite_quadratic_module
- psage.modform.weilrep_tools

I will try to keep this repository in sync with my psage fork.
You can manually sync with my psage fork by applying the following
commands in the repository's top level directory:

```
  git remote add sehlen-psage https://github.com/sehlen/psage
  git subtree pull -p psage_subtree sehlen-psage devel
```
