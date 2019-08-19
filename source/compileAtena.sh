#! /bin/sh

# pwd
# >> hoge/source
# ./compile.sh cmake
# >> try compile with cmake

if [ "$TestArea" = "" ] ; then 
  echo "TestArea is missing"
  exit 1
fi
echo "TestArea = "$TestArea

BUILD_DIR=$TestArea/../build
rm -r $BUILD_DIR
if [ -d "$BUILD_DIR" ] ; then
  echo "Directory for build is $BUILD_DIR"
else
  echo "Creating $BUILD_DIR ..." && mkdir -p $BUILD_DIR
fi

cd $BUILD_DIR
if [ $# -gt 0 ] && [ $1 = "cmake" ]; then
  cmake $TestArea/athena/Projects/WorkDir
fi
make -j8 -C $BUILD_DIR
cd -
