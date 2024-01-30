## compilation
```sh
git clone --recurse-submodules https://github.com/clairekth/Snowman.git
cd Snowman
git submodule update --init
mkdir build
cd build
cmake ..  
make
./snowman
```
