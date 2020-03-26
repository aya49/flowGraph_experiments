/*
 * Implemented by Adrin Jalali <adrin.jalali@gmail.com>
 * October 2011
 * Licensed under GPLv3
*/
#include "numberGenerator.h"
#include <cstring>
#include <cstdio>

#include "exceptions.h"

using namespace tfl;

NumberGenerator::NumberGenerator(int len, int base, int k): len(len), base(base), k(k), started(false)
{
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	for (int i = 0; i < len; i++)
	{
		this->v_base[i] = base;
	}
}

NumberGenerator::NumberGenerator(int len, int base, const char *text): len(len), base(base), started(true)
{
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	for (int i = 0; i < len; i++)
	{
		this->_data[i] = text[i] - '0';
		this->v_base[i] = base;
	}
	this->k = nonZeroCount();
}

NumberGenerator::NumberGenerator(int len, const int *base, int k): len(len), k(k), started(false)
{
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	memcpy(this->v_base, base, sizeof(int) * len);
}

NumberGenerator::NumberGenerator(int len, const int *base, const char *text): len(len), started(true)
{
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	memcpy(this->v_base, base, sizeof(int) * len);
	for (int i = 0; i < len; i++)
	{
		this->_data[i] = text[i] - '0';
	}
	this->k = nonZeroCount();
}

NumberGenerator::NumberGenerator(const NumberGenerator &other)
{
	len = other.len;
	base = other.base;
	k = other.k;
	started = other.started;
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	memcpy(this->_data, other._data, len);
	memcpy(this->_text, other._text, len + 1);
	memcpy(this->v_base, other.v_base, sizeof(int) * len);
}

NumberGenerator &NumberGenerator::operator=(const NumberGenerator &other)
{
	len = other.len;
	base = other.base;
	k = other.k;
	started = other.started;
	this->_data = new char[len];
	this->_text = new char[len + 1];
	this->v_base = new int[len];
	memcpy(this->_data, other._data, len);
	memcpy(this->_text, other._text, len + 1);
	memcpy(this->v_base, other.v_base, sizeof(int) * len);
	return *this;
}

NumberGenerator::~NumberGenerator()
{
	delete[] _data;
	delete[] _text;
	delete[] v_base;
}

const char *NumberGenerator::data()
{
	return this->_data;
}

const char *NumberGenerator::text()
{
	dataToText(_data, _text);
	return _text;
}

/*****************************************************************************
 * Get first phenotype in sequence 
 *
 ****************************************************************************/


NumberGenerator &NumberGenerator::first()
{
	memset(_data, 0, len);
	for (int i = 0; i < k; i++) 
		_data[i] = 1;

	return *this;
}

/*****************************************************************************
 * Get next phenotype in sequence 
 *
 ****************************************************************************/


NumberGenerator &NumberGenerator::next()
{
	if (!started)
	{
		started = true;
		return first();
	}

	int inc = firstIncreasableIndex();
	//Rprintf("%d !!!\n", inc);
	if (inc != -1)
	{
		if (inc > 0 && !this->_data[inc])
		{
			//Rprintf("happened! inc:%d\n", inc);
			this->_data[inc - 1] = 0;
		}
		this->_data[inc]++;
		resetLessSignificantNonZeroDigits(inc);
	}
	else
	{
		int zinc = firstZeroGreaterThanANonZeroIndex();
		if (zinc == -1)
			throw Exception("no more to generate!");

		this->_data[zinc] = 1;
		if (zinc > 0)
			this->_data[zinc - 1] = 0;
		resetLessSignificantNonZeroDigits(zinc);
	}
	return *this;
}

/*****************************************************************************
 * Are we at the last phenotype? 
 *
 ****************************************************************************/

bool NumberGenerator::hasNext()
{
	if (!started)
		return true;
	if (firstIncreasableIndex() == -1 && firstZeroGreaterThanANonZeroIndex() == -1)
		return false;
	return true;
}

/*****************************************************************************
 * 
 *
 ****************************************************************************/

void NumberGenerator::resetLessSignificantNonZeroDigits(int index)
{
	int count = 0;
	for (int i = 0; i < index; i++)
		if (this->_data[i])
			count++;

	//Rprintf("count: %d ", count);
	for (int i = 0; i < count; i++)
	{
		//Rprintf("%d ", i);
		this->_data[i] = 1;
	}
	for (int i = count; i < index; i++)
	{
		//Rprintf("_%d ", i);
		this->_data[i] = 0;
	}
	//Rprintf("\n");
}

int NumberGenerator::greatestNonZeroIndex()
{
	for (int i = len - 1; i >= 0; i--)
		if (this->_data[i])
			return i;
	return -1;
}

int NumberGenerator::firstIncreasableIndex()
{
	bool non_zero_seen = false;
	for (int i = 0; i < len; i++)
	{
		if (this->_data[i])
			non_zero_seen = true;
		if (this->_data[i] < v_base[i] - 1 && (this->_data[i] || non_zero_seen))
			return i;
	}
	return -1;
}

int NumberGenerator::firstZeroGreaterThanANonZeroIndex()
{
	for (int i = 0; i < this->len - 1; i++)
		if (this->_data[i] && !this->_data[i + 1])
			return i + 1;
	return -1;
}

int NumberGenerator::getLength()
{
	return len;
}


/*****************************************************************************
 * 
 *
 ****************************************************************************/

NumberGenerator NumberGenerator::neighbor(int n_k)
{
	NumberGenerator result(len, v_base, k - 1);
	memcpy(result._data, _data, len);
	int c = 0;
	for (int i = 0; i < len; i++)
	{
		if (result._data[i])
			c++;
		if (c - 1 == n_k)
		{
			result._data[i] = 0;
			return result;
		}
	}
	char msg[1000];
	sprintf(msg, "no neighbor at non zero digit given. number: %s\t input: %d", this->text(), n_k);
	throw Exception(msg);
}

void NumberGenerator::dataToText(const char *data, char *text)
{
	for (int i = 0; i < len; i++)
		text[i] = data[i] + '0';
	text[len] = 0;
}

int NumberGenerator::nonZeroCount()
{
	int c = 0;
	for (int i = 0; i < len; i++)
		if (_data[i])
			c++;
	return c;
}

int NumberGenerator::firstDiff(const NumberGenerator &other)
{
	for (int i = 0; i < len; i++)
		if (this->_data[i] != other._data[i])
			return i;
	return -1;
}
