/*
 * Implemented by Adrin Jalali <adrin.jalali@gmail.com>
 * October 2011
 * Licensed under GPLv3
*/
#ifndef _NUMBER_GENERATOR_
#define _NUMBER_GENERATOR_

namespace tfl
{

	class NumberGenerator
	{
		private:
			char *_data;
			char *_text;
			int len, base, k;
			int* v_base;
			bool started;

			int firstIncreasableIndex();
			int firstZeroGreaterThanANonZeroIndex();
			//reset all nonzero less significant digits than "index" to 1
			void resetLessSignificantNonZeroDigits(int index);
			void dataToText(const char *data, char *text);

			NumberGenerator() { }

		public:
			NumberGenerator(const NumberGenerator &other);
			NumberGenerator &operator=(const NumberGenerator &other);

			NumberGenerator(int len, int base, int k);
			NumberGenerator(int len, int base, const char *text);
			NumberGenerator(int len, const int *base, int k);
			NumberGenerator(int len, const int *base, const char *text);

			NumberGenerator &first();
			NumberGenerator &next();
			bool hasNext();
			int nonZeroCount();
			int getLength();

			int greatestNonZeroIndex();
			int firstDiff(const NumberGenerator &other);

			const char *text();
			const char *data();
			//zero based neighbor number
			NumberGenerator neighbor(int n_k);
			
			~NumberGenerator();
	};
}


#endif
