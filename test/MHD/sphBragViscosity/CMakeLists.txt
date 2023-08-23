enable_idefix_property(Idefix_MHD)
## The current test is not physics, but merely a non-regression check
## The physical test needs to have the divergence of the velocity field artifiacially set to zero
#replace_idefix_source(fluid/braginskii/bragViscosity.cpp noDivBragViscosity.cpp)
