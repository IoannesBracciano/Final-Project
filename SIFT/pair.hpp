#pragma once
/** The next portion of code was ported from the matlab interface of the vlfeat library
as there was no suitable code for matching between two image descriptors in C API
and was copied and pasted almost "as is".
See http://stackoverflow.com/questions/26925038/how-to-use-vlfeat-sift-matching-function-in-c-code
*/

/// <summary>Holds matched descriptor pairs</summary>
typedef struct {
	/// <summary>Descriptor of first image</summary>
	int k1;
	/// <summary>Descriptor of second image</summary>
	int k2;
	/// <summary>Score of matching</summary>
	double score;
} Pair;

/// <summary>Compares the descriptors of two different images and matches the pairs</summary>
/// <param name="pairs">pointer to a vector to be filled with the matched pairs.
///                     Needs to be big enough to hold all possible pairs</param>
/// <param name="descr1">pointer to the descriptor matrix of the first image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="descr2">pointer to the descriptor matrix of the second image,
///                      in row major order (K1 rows x ND columns)</param>
/// <param name="K1">Number of descriptors of first image</param>
/// <param name="K2">Number of descriptors of second image</param>
/// <param name="ND">Descriptor size (=128 for SIFT)</param>
/// <param name="thresh">Ratio test threshold value (=1.5)</param>
/// <returns>a pointer pointing to the end of the vector conaitning the matched pairs</returns>
Pair * compare(Pair * pairs, const float *descr1, const float *descr2,
	int K1, int K2, int ND = 128, float thresh = 1.5)
{
	Pair * _pairs = pairs; // Copy pointer to avoid changing the original
	int k1, k2;

	/* Loop over 1st image descr. */
	for (k1 = 0; k1 < K1; ++k1, descr1 += ND) {
		float best = FLT_MAX;
		float second_best = FLT_MAX;
		int bestk = -1;

		/* Loop over 2nd image descr. and find the 1st and 2nd closest descr. */
		for (k2 = 0; k2 < K2; ++k2, descr2 += ND) {
			int bin;
			float acc = 0;

			/* Compute the square L2 distance between descriptors */
			for (bin = 0; bin < ND; ++bin) {
				float delta = descr1[bin] - descr2[bin];
				acc += delta*delta;
				if (acc >= second_best)
					break;
			}

			if (acc < best) {
				second_best = best;
				best = acc;
				bestk = k2;
			}
			else if (acc < second_best) {
				second_best = acc;
			}
		}

		/* Rewind */
		descr2 -= ND*K2;

		/* Record the correspondence if the best descr. passes the ratio test */
		if (thresh * best < second_best && bestk != -1) {
			_pairs->k1 = k1;
			_pairs->k2 = bestk;
			_pairs->score = best;
			_pairs++;
		}
	}

	return _pairs;
}
/** End of ported code section */