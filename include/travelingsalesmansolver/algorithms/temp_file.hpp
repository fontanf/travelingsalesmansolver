#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <cstdio>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#endif

namespace travelingsalesmansolver
{

/**
 * Reserve a unique temporary file path.
 *
 * On Linux, this prefers a tmpfs ('/dev/shm') so that the instance /
 * parameter / solution files exchanged with external solver subprocesses
 * (LKH, Concorde, ...) stay in memory instead of hitting disk, falling back
 * to '$TMPDIR' / '/tmp' otherwise (e.g. macOS, which has no '/dev/shm').
 * On Windows, it uses the regular temporary directory ('GetTempPath'); there
 * is no standard memory-backed equivalent available there.
 *
 * The file is created (to guarantee uniqueness, unlike 'tmpnam') and then
 * removed immediately, so the caller gets a name that is free to use,
 * consistently with how callers create the actual files (either via
 * 'ofstream' or by letting the external subprocess write them).
 */
inline std::string make_temp_path(const std::string& prefix)
{
#if defined(_WIN32)
    char dir[MAX_PATH];
    if (GetTempPathA(MAX_PATH, dir) == 0) {
        throw std::runtime_error("Unable to determine a temporary directory.");
    }

    // 'GetTempFileNameA' only uses the first three characters of 'prefix'.
    char path[MAX_PATH];
    if (GetTempFileNameA(dir, prefix.c_str(), 0, path) == 0) {
        throw std::runtime_error(
                "Unable to create a temporary file in \"" + std::string(dir) + "\".");
    }
    std::remove(path);
    return std::string(path);
#else
    struct stat info;
    std::string dir = "/dev/shm";
    if (stat(dir.c_str(), &info) != 0 || !S_ISDIR(info.st_mode)) {
        const char* tmpdir = std::getenv("TMPDIR");
        dir = (tmpdir != nullptr) ? tmpdir : "/tmp";
    }

    std::string path_template = dir + "/" + prefix + "XXXXXX";
    std::vector<char> buffer(path_template.begin(), path_template.end());
    buffer.push_back('\0');

    int fd = mkstemp(buffer.data());
    if (fd == -1) {
        throw std::runtime_error(
                "Unable to create a temporary file in \"" + dir + "\".");
    }
    close(fd);
    std::string path(buffer.data());
    std::remove(path.c_str());
    return path;
#endif
}

}
