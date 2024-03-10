load("//travelingsalesmansolver:variables.bzl", "STDCPP")

cc_library(
        name = "distances",
        hdrs = [
                "commons.hpp",
                "distances.hpp",
                "distances_builder.hpp",
                "distances_explicit.hpp",
                "distances_explicit_triangle.hpp",
                "distances_euc_2d.hpp",
                "distances_ceil_2d.hpp",
                "distances_geo.hpp",
                "distances_att.hpp",
        ],
        srcs = [
                "distances.cpp",
                "distances_builder.cpp",
        ],
        copts = STDCPP,
        deps = [
                "@optimizationtools//optimizationtools/utils:utils",
        ],
        visibility = ["//visibility:public"],
)
