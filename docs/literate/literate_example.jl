# # Code Example
#
# This is an example of adding a code example compiled with Literate.jl in the docs.
#
# First we want to load our package with `using`

using CellularSheaves
using LinearAlgebra

# ## Using the Package
# Here is a simple example of creating a sheaf using the `@cellular_sheaf` macro.

A = [1.0 0.0 1.0 0.0]
B = [1.0 0.0 0.0 1.0]
C = [1.0 0.0 0.0 0.0]

sheaf = @cellular_sheaf A, B, C begin
    x::Stalk{4}, y::Stalk{4}, z::Stalk{4}

    A(x) == B(y)
    A(x) == C(z)
    B(y) == C(z)

end

# Let's compute a global section of this sheaf.
# We start with a random 0-cochain:
x0 = rand(sum(vertex_stalks(sheaf)))

global_section = nearest_global_section(sheaf, x0)

# Now we can check that this is indeed a global section by
# verifying that the coboundary map applied to it is zero.
d = coboundary_map(sheaf)

norm(d * global_section)
