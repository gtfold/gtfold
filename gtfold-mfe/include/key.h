class Key {
	public:
	Key(unsigned char type, int i, int j);
	unsigned char type;
	int i;
	int j;
	unsigned char getType();
	int getI();
	int getJ();
	bool operator<(const Key& other) const;
	bool operator==(const Key& other) const;
	
};
