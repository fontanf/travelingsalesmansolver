load("//travelingsalesmansolver:variables.bzl", "STDCPP")

cc_library(
        name = "distances",
        hdrs = [
                "distances.hpp",
                "distances_builder.hpp",
        ],
        srcs = [
                "distances.cpp",
                "distances_builder.cpp",
        ],
        deps = [
                "@optimizationtools//optimizationtools/utils:utils",
        ],
        visibility = ["//visibility:public"],
)


cc_library(
        name = "travelingsalesmansolver",
        hdrs = [
                "instance.hpp",
                "solution.hpp",
                "algorithm_formatter.hpp",
        ],
        srcs = [
                "instance.cpp",
                "solution.cpp",
                "algorithm_formatter.cpp",
        ],
        deps = [
                ":distances",
                "@optimizationtools//optimizationtools/utils:utils",
                "@optimizationtools//optimizationtools/containers:containers",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "main",
        srcs = ["main.cpp"],
        deps = [
                "//travelingsalesmansolver/algorithms:lkh",
                "//travelingsalesmansolver/algorithms:concorde",
                "@boost//:program_options",
        ],
        linkopts = select({
                    "@bazel_tools//src/conditions:windows": [],
                    "//conditions:default":                 ["-lpthread"],
        }),
)

