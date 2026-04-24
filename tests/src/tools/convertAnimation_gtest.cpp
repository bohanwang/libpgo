#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef PGO_TEST_CONVERT_ANIMATION_BIN
#define PGO_TEST_CONVERT_ANIMATION_BIN ""
#endif

namespace
{
namespace fs = std::filesystem;

std::string quotePath(const fs::path &path)
{
  return "\"" + path.string() + "\"";
}

std::string shellExecutable(const fs::path &path)
{
#ifdef _WIN32
  return "call " + quotePath(path);
#else
  return quotePath(path);
#endif
}

fs::path getConvertAnimationBinaryPath()
{
  const char *envBin = std::getenv("CONVERT_ANIMATION_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return fs::path(envBin);

  if (std::string(PGO_TEST_CONVERT_ANIMATION_BIN).empty())
    return {};

  return fs::path(PGO_TEST_CONVERT_ANIMATION_BIN);
}

int runCommand(const std::string &command)
{
  return std::system(command.c_str());
}

class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto base = fs::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = base / ("libpgo-convertAnimation-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

      std::error_code ec;
      if (fs::create_directories(path_, ec))
        return;
    }

    throw std::runtime_error("Failed to create a unique temporary directory");
  }

  ~ScopedTempDir()
  {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }

  const fs::path &path() const { return path_; }

private:
  fs::path path_;
};

void writeTextFile(const fs::path &path, const std::string &contents)
{
  fs::create_directories(path.parent_path());
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());
  out << contents;
}
}  // namespace

TEST(ConvertAnimationCli, DefaultsOutputFolderToConfigDirectory)
{
  const fs::path binary = getConvertAnimationBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configDir = tempDir.path() / "anim";
  const fs::path configPath = configDir / "anim.json";
  const fs::path meshPath = configDir / "mesh.obj";
  const fs::path frame0Path = configDir / "frames" / "frame0000.obj";
  const fs::path frame1Path = configDir / "frames" / "frame0001.obj";
  const fs::path expectedAbcPath = configDir / "tri.abc";

  writeTextFile(meshPath, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame0Path, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame1Path, "v 0 0 0\nv 1.2 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(configPath, R"({
  "save-cache": 0,
  "meshes": [
    {
      "name": "tri",
      "driving-mesh": "mesh.obj",
      "sequence": "frames/frame{:04d}.obj",
      "sequence-type": "objmesh",
      "sequence-range": [0, 2]
    }
  ]
})");

  std::ostringstream command;
  command << shellExecutable(binary)
          << " " << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(expectedAbcPath));
  ASSERT_GT(fs::file_size(expectedAbcPath), 0);
}

TEST(ConvertAnimationCli, RespectsExplicitOutputFolder)
{
  const fs::path binary = getConvertAnimationBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configDir = tempDir.path() / "anim";
  const fs::path outputDir = tempDir.path() / "abc-out";
  const fs::path configPath = configDir / "anim.json";
  const fs::path meshPath = configDir / "mesh.obj";
  const fs::path frame0Path = configDir / "frames" / "frame0000.obj";
  const fs::path frame1Path = configDir / "frames" / "frame0001.obj";
  const fs::path expectedAbcPath = outputDir / "tri.abc";

  fs::create_directories(outputDir);
  writeTextFile(meshPath, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame0Path, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame1Path, "v 0 0 0\nv 1.2 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(configPath, R"({
  "save-cache": 0,
  "meshes": [
    {
      "name": "tri",
      "driving-mesh": "mesh.obj",
      "sequence": "frames/frame{:04d}.obj",
      "sequence-type": "objmesh",
      "sequence-range": [0, 2]
    }
  ]
})");

  std::ostringstream command;
  command << shellExecutable(binary)
          << " " << quotePath(configPath)
          << " " << quotePath(outputDir);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(expectedAbcPath));
  ASSERT_GT(fs::file_size(expectedAbcPath), 0);
}
