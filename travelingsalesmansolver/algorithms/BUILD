cc_library(
        name = "concorde",
        hdrs = ["concorde.hpp"],
        srcs = ["concorde.cpp"],
        deps = ["//travelingsalesmansolver:travelingsalesmansolver"],
        visibility = ["//visibility:public"],
)

cc_library(
        name = "lkh",
        hdrs = ["lkh.hpp"],
        srcs = ["lkh.cpp"],
        deps = ["//travelingsalesmansolver:travelingsalesmansolver"],
        visibility = ["//visibility:public"],
)

cc_library(
        name = "algorithms",
        hdrs = ["algorithms.hpp"],
        srcs = ["algorithms.cpp"],
        deps = [
                ":concorde",
                ":lkh",
                "@boost//:program_options",
        ],
        visibility = ["//visibility:public"],
)

