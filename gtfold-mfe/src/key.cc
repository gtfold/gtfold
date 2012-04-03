#include "key.h"

unsigned char type;
int i;
int j;

Key::Key(unsigned char type, int i, int j) {
	this->type = type;
	this->i = i;
	this->j = j;
}

bool Key::operator<(const Key& other) const {
	return (this->type != other.type ? this->type < other.type : this->i != other.i ? this->i < other.i : this->j < other.j);
}

bool Key::operator==(const Key& other) const {
	return (this->type == other.type && this->i == other.i && this->j == other.j);
}

unsigned char Key::getType() {
	return type;
}

int Key::getI() {
	return i;
}

int Key::getJ() {
	return j;
}
