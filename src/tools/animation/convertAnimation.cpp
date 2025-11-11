#include "animationLoader.h"
#include "pgoLogging.h"
#include "initPredicates.h"

int main(int argc, char *argv[])
{
  pgo::Mesh::initPredicates();
  pgo::Logging::init();

  pgo::AnimationIO::AnimationLoader loader;
  if (loader.load(argv[1]) != 0) {
    return 1;
  }

  loader.saveABC(argv[2]);

  return 0;
}