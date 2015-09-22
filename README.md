# JLDArchives

[JLD](https://github.com/JuliaLang/JLD.jl) is an HDF5-based file format for storing data for the Julia language.

This is a repository of old *.jld files, useful for ensuring that we preserve backwards compatibility.

At present, there is no "runtime" code in this repository: it merely includes *.jld files and some test scripts.
At some point if it becomes difficult to maintain backwards compatibility, this might become a home
for converting old formats to more modern ones.
(For example, by archiving old versions of the JLD format implementation, and then saving all the same
variables using the modern format.)

To test backwards compatility, simply say `Pkg.test("JLDArchives")` from julia,
or run the `runtests.jl` script from inside the `test` folder.
