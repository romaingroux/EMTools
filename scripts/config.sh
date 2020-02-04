
dir=$(pwd)
sed -i "s@<path_to_project_root>@$dir@g"          build.sh
sed -i "s@<path_to_project_root>@$dir@g"     CMakeLists.txt
sed -i "s@<path_to_project_root>@$dir@g" src/CMakeLists.txt
