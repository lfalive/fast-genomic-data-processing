#include "bitmap.h"
#include <cassert>
#include <iostream>

bitmap::bitmap(uint64_t _size): size(_size){
  uint64_t lines = (size + 63) / 64;
  array = new std::atomic_uint64_t[lines];
  for(uint64_t i=0; i < lines; i++){
    array[i] = 0;
  }
  std::cout << "bitmap size: " << double(lines) * sizeof(uint64_t) / 1024 / 1024 / 1024 <<"GB" <<std::endl;
}

bitmap::~bitmap(){
  delete[] array;
}

void bitmap::set(uint64_t pos){
  if (pos >= size) {
  	std::cerr << "set Error: Position " << pos << " out of bounds (size: " << size << ")" << std::endl;
        assert(pos < size);
  }
  uint64_t line = pos >> 6; // line = pos / 64
  uint64_t offset = pos % 64;
  uint64_t value_or = 1ll << offset;
  array[line] |= value_or;
}

bool bitmap::get(uint64_t pos){
  if (pos >= size) {
        std::cerr << "get Error: Position " << pos << " out of bounds (size: " << size << ")" << std::endl;
        assert(pos < size);
  }
  uint64_t line = pos >> 6;
  uint64_t offset = pos % 64;
  uint64_t line_value = array[line];
  uint64_t value_or = 1ll << offset;
  return (value_or & line_value) != 0;
}
