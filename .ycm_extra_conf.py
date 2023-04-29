def Settings(**kwargs):
    return {
            'flags': [
                '-x', 'c++',
                '-Wall', '-Wextra', '-Werror',
                '-I', '.',

                '-I', './bazel-travelingsalesmansolver/external/'
                'json/single_include/',

                '-I', './bazel-travelingsalesmansolver/external/'
                'googletest/googletest/include/',

                '-I', './bazel-travelingsalesmansolver/external/'
                'boost/',

                # optimizationtools
                '-I', './bazel-travelingsalesmansolver/external/'
                # '-I', './../'
                'optimizationtools/',

                # localsearchsolver
                '-I', './bazel-travelingsalesmansolver/external/'
                # '-I', './../'
                'localsearchsolver/',

                # treesearchsolver
                '-I', './bazel-travelingsalesmansolver/external/'
                # '-I', './../'
                'treesearchsolver/',
                ],
            }
