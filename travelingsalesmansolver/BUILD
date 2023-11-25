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
        ],
        srcs = [
                "instance.cpp",
                "solution.cpp",
        ],
        deps = [
                ":distances",
                "@optimizationtools//optimizationtools/utils:info",
                "@optimizationtools//optimizationtools/utils:utils",
                "@optimizationtools//optimizationtools/containers:containers",
        ],
        visibility = ["//visibility:public"],
)

cc_binary(
        name = "main",
        srcs = ["main.cpp"],
        deps = [
                "//travelingsalesmansolver/algorithms:algorithms",
        ],
        linkopts = select({
                    "@bazel_tools//src/conditions:windows": [],
                    "//conditions:default":                 ["-lpthread"],
        }),
)

