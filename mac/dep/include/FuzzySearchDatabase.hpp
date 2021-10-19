#pragma once

#ifndef FUZZY_SEARCH_DATABASE_HPP
#define FUZZY_SEARCH_DATABASE_HPP


// Copyright 2020 Nils Jonas Norberg <jnorberg@gmail.com>
// license: BSD-3-Clause ( https://opensource.org/licenses/BSD-3-Clause )

// updates, requests, comments:
// https://bitbucket.org/j_norberg/fuzzysearchdatabase



#include <vector>
#include <unordered_map>
#include <string>
#include <cctype>
#include <cstring>
#include <algorithm>

namespace fuzzysearch
{

	// helper class for the fuzzy distance
	class HelperFunctions
	{
		static inline int min3(int a, int b, int c)
		{
			return (a < b) ? (a < c ? a : c) : (b < c ? b : c);
		}

		static inline int max(int a, int b)
		{
			return a > b ? a : b;
		}

		static inline float max(float a, float b)
		{
			return a > b ? a : b;
		}

		// 
		static inline float substringScore(const char* qw, int qwL, const char* w, int wL)
		{
			int subStrEnd = wL - qwL;
			// maybe limit the lengths here? 

			for (int substrPos = 0; substrPos <= subStrEnd; ++substrPos)
			{
				if (0 == memcmp(w + substrPos, qw, qwL))
				{
					// score substrings close to the beginning of the word higher
					// score higher if length of substring is closer to the whole word
					float looseFit = (float)(wL - qwL);
					float score = 20.0f / (20.0f + substrPos + looseFit * 0.5f);
					return score;
				}
			}

			return 0.0f;
		}

		static inline int levDistance(const char* a, int aL, const char* b, int bL)
		{
			// ------------------------------------------------------
			// https://en.wikipedia.org/wiki/Levenshtein_distance
				
			// common prefix
			while (aL > 0 && bL > 0 && a[0] == b[0])
			{
				++a;
				++b;
				--aL;
				--bL;
			}

			// common suffix
			while (aL > 0 && bL > 0 && a[aL - 1] == b[bL - 1])
			{
				--aL;
				--bL;
			}

			// simple case
			if (aL < 1) return bL;
			if (bL < 1) return aL;

			enum { kBuf = 16, kLen = kBuf - 1 };

			// clamp length
			if (aL > kLen) aL = kLen;
			if (bL > kLen) bL = kLen;

			// used during distance-calculations
			// currently have a length-limit on words, should be ok
			// this could in theory be just two rows ( or even one )
			int _dist[kBuf * kBuf];

			int stride = aL + 1;

			// init top
			for (int x = 0; x <= aL; x++)
				_dist[x] = x;

			// fill rest of matrix
			for (int y = 1; y <= bL; y++)
			{
				int index = y * stride + 1;

				// init-left side
				_dist[index - 1] = y;

				int bChar = b[y - 1];

				// in theory, could unroll this a little-bit
				for (int x = 1; x <= aL; ++x, ++index)
				{
					int substitutionCost = (a[x - 1] == bChar) ? 0 : 1;

					int t = min3(
						_dist[index - 1] + 1,
						_dist[index - stride] + 1,
						_dist[index - stride - 1] + substitutionCost
					);

					_dist[index] = t; // write
				}
			}

			// return
			int dist = _dist[bL * stride + aL];
			return dist;
		}

		static inline bool isDivider(int c)
		{
			// possible improvement, support utf8
			// maybe this is a bit agressive
			// this will treat every character not mentioned below as a delimiter
			if ('a' <= c && c <= 'z') return false;
			if ('A' <= c && c <= 'Z') return false;
			if ('0' <= c && c <= '9') return false;
			return true;
		}



	public:

		static inline void toLower( std::string& s )
		{
			for ( size_t i = 0 ; i < s.size() ; ++i )
				s[i] = (char)std::tolower( s[i] );
		}

		static inline std::vector< std::string > splitString( std::string s )
		{
			std::vector< std::string > ss;

			size_t left = 0;
			size_t right = 0;

			for (; left < s.size(); ++left)
			{
				// 1. skip all dividers
				if (isDivider(s[left]))
					continue;

				right = left + 1;

				// 2. skip all non-dividers 
				for (; right < s.size(); ++right)
					if (isDivider(s[right]))
						break;

				ss.push_back(s.substr(left, right - left));
				left = right;
			}

			return ss;
		}



		// qw - query-Word
		// this gives a score on how well a query-word fits a word
		// 1.0 is a full match
		static float scoreQueryWordToWord(const char* qw, int qwL, const char* w, int wL)
		{
			float subStrScore = 0.0f;

			if (qwL == wL)
			{
				// test for full match
				if (0 == memcmp(qw, w, qwL))
					return 1.0f;
			}
			else if (qwL < wL)
			{
				// this could be a substring if the query is shorter or same
				subStrScore = substringScore(qw, qwL, w, wL);
			}

			// skip fuzzy-calculation if word is "much" longer than query-word
			const int kLongerLim = 4;
			if (wL > qwL + kLongerLim)
				return subStrScore;

			// useless to do fuzzy on single character
			if (qwL < 2)
				return subStrScore;

			// if needed use fuzzy
			// could build in the distance-limit here to improve early-out
			int fuzzyDist = levDistance(qw, qwL, w, wL);

			// if the distance is very high ( more than half of the query-word-length )
			// it's not good enough
			int distanceLim = (qwL + 1) / 2;
			if (fuzzyDist >= distanceLim)
				return subStrScore;

			float fuzzyDistFrac = (float)fuzzyDist / (float)distanceLim;
			float fuzzyScore = 1.0f - fuzzyDistFrac;

			float finalScore = max(fuzzyScore, subStrScore);
			return finalScore;
		}


	};










	template <typename Key = std::string>
	class Database
	{
	public:

		// ------------------------------------------------------
		// public struct, used for returning the results

		struct Result
		{
			Key key;
			float score;
		};

	private:


		// helper struct to store words close in memory
		struct WordStorage
		{
			std::vector<char> _wordData; // all characters from a word, in order in memory
			std::vector<int> _wordEnd; // the "end" marker for each word
			std::unordered_map< std::string, int > _wordMap;

			void clear()
			{
				_wordData.clear();
				_wordEnd.clear();
				_wordMap.clear();
			}

			int addWord(std::string word)
			{
				// if exist, return index
				auto found = _wordMap.find(word);
				if (found != _wordMap.end())
				{
					return found->second;
				}

				// add to data
				_wordData.insert(_wordData.end(), word.begin(), word.end());

				// add end-marker
				int wordIndex = (int)_wordEnd.size();
				_wordEnd.push_back((int)_wordData.size());

				// add to map
				_wordMap.insert({ word, wordIndex });

				// return word-index
				return wordIndex;
			}
		};

		//
		struct WordFromField {
			size_t _wordIndex;
			size_t _fieldIndex;

			WordFromField(size_t wIndex, size_t fIndex)
				: _wordIndex(wIndex)
				, _fieldIndex(fIndex)
			{
			}

			bool operator<(const WordFromField& other) const {
				if (_wordIndex != other._wordIndex)
					return _wordIndex < other._wordIndex;

				return _fieldIndex < other._fieldIndex;
			}

			bool operator==(const WordFromField& other) const {
				return
					_wordIndex == other._wordIndex &&
					_fieldIndex == other._fieldIndex;
			}
		};

		// temporary struct to make sorting smooth
		struct TempResult
		{
			size_t _entryIndex;
			float _score;

			TempResult()
			{
				set(0.0f, 0);
			}

			TempResult(float score, size_t entryIndex)
			{
				set(score, entryIndex);
			}

			void set(float score, size_t entryIndex)
			{
				_score = score;
				_entryIndex = entryIndex;
			}

			bool operator<(const TempResult& other) const {
				return _score > other._score; // reversed to sort high score on top
			}
		};

		// this is the entry, all words are stored in the "WordStorage", only indices here
		struct Entry
		{
			Key _key;
			std::vector<WordFromField> _words;	// index into words from strings
		};


		// ------------------------------------------------------

		// memory-friendly storage for all words
		WordStorage _ws;

		// all entries
		std::vector<Entry> _entries;

		//
		std::vector<float> _fieldWeights;

		//
		float _threshold;

		// "temporary" used during search
		mutable std::vector<float> _scorePerWord;
		mutable std::vector< TempResult > _tempResults;


		static void scoreEveryWord(std::vector<float>& scores, const WordStorage& ws, const std::string& queryWord)
		{
			scores.resize(ws._wordEnd.size());

			const char* qw = queryWord.c_str();
			const int qwL = (int)queryWord.size();

			// skip words much shorter than query-word
			int lim = qwL < 4 ? qwL - 1 : qwL - 2;

			// visit every word in memory-order
			int bI = 0;
			for (size_t i = 0; i < ws._wordEnd.size(); ++i)
			{
				int eI = ws._wordEnd[i];
				const char* b = &ws._wordData[bI];
				int bL = eI - bI;

				// skip short strings when query is longer
				if (bL < lim)
					scores[i] = 0.0f;
				else
					scores[i] = HelperFunctions::scoreQueryWordToWord(qw, qwL, b, bL);

				bI = eI;
			}
		}

		float scoreEntry(const Entry& e) const
		{
			float score = 0;

			// for each word in this entry
			// find highest match
			for (size_t ei = 0; ei < e._words.size(); ++ei)
			{
				const WordFromField& wff = e._words[ei];
				float localScore = _scorePerWord[wff._wordIndex];
				float weight = _fieldWeights[wff._fieldIndex];
				localScore *= weight;

				if (localScore > score)
					score = localScore;
			}

			return score;
		}

		void scoreEveryEntry(bool first) const
		{
			if (first)
			{
				// set
				for (size_t i = 0; i < _tempResults.size(); ++i)
					_tempResults[i].set(scoreEntry(_entries[i]), i);
			}
			else
			{
				// multiply
				for (size_t i = 0; i < _tempResults.size(); ++i)
					_tempResults[i]._score *= scoreEntry(_entries[i]);
			}
		}

	public:


		Database()
		{
			reset();
		}

		void reset()
		{
			_ws.clear();
			_scorePerWord.clear();
			_tempResults.clear();
			_entries.clear();
			_fieldWeights.clear();
			_threshold = 0.1f;
		}

		void addEntry( Key key, const std::vector<std::string>& fields )
		{
			// ensure we have enough weights
			while ( _fieldWeights.size() < fields.size() )
				_fieldWeights.push_back( 1.0f );

			// create entry
			Entry e;
			e._key = key;
			std::vector<WordFromField>& eWs = e._words;

			// iterate fields
			for ( size_t fieldIndex = 0; fieldIndex < fields.size(); ++fieldIndex )
			{
				auto words = HelperFunctions::splitString(fields[fieldIndex]);

				// iterate words in this field
				for (size_t wi = 0; wi < words.size(); ++wi)
				{
					std::string& word = words[wi];
					HelperFunctions::toLower(word);

					// add the word ( de-dup happens here )
					int wordIndex = _ws.addWord(word);
					eWs.emplace_back( wordIndex, fieldIndex );
				}
			}

			// key will never be found if no strings, so don't add
			if (eWs.size() < 1)
				return;

			// sort indices ( in memory-order )
			std::sort(eWs.begin(), eWs.end());

			// we remove duplicates, ( if we knew about the weights here, we could remove duplicates across fields too... )
			eWs.erase(std::unique(eWs.begin(), eWs.end()), eWs.end());

			// finally add
			_entries.push_back(e);
		}

		// each field can have a weight (defaults to 1)
		void setWeights( const std::vector<float>& fieldWeights )
		{
			// ensure we have enough weights
			while (_fieldWeights.size() < fieldWeights.size())
				_fieldWeights.push_back(1.0f);

			for (size_t i = 0; i < fieldWeights.size(); ++i)
				_fieldWeights[i] = fieldWeights[i];
		}

		// any search-result scoring below this will not be returned from the search
		void setThreshold(float threshold)
		{
			_threshold = threshold;
		}

		// returns a big list of results ( sorted and with a score )
		std::vector< Result > search( std::string queryString ) const
		{
			std::vector< Result > resultsWithKey;

			// mirrors all entries
			_tempResults.resize(_entries.size());

			// 0. prepare query (query-string -> query-words)
			auto qws = HelperFunctions::splitString(queryString);
			for (size_t i = 0; i < qws.size(); ++i)
				HelperFunctions::toLower(qws[i]);

			// 1. loop over each word in query
			for (size_t qi = 0; qi < qws.size(); ++qi)
			{
				const std::string& qWord = qws[qi];

				// 2. score every word against this query-word
				scoreEveryWord(_scorePerWord, _ws, qWord);

				// 3. score each entry
				scoreEveryEntry(qi == 0);
			}

			// at this point all scores are in tempResults vector
			// only sorting left

			// erase below threshold
			_tempResults.erase(std::remove_if(_tempResults.begin(), _tempResults.end(), [&](TempResult tr) { return tr._score < _threshold; }), _tempResults.end());

			// sort all that remain
			std::sort( _tempResults.begin(), _tempResults.end() );

			// finally copy to the result vector
			resultsWithKey.resize( _tempResults.size() );
			for ( size_t i = 0; i < _tempResults.size(); ++i )
			{
				const TempResult& src = _tempResults[i];
				Result& dst = resultsWithKey[i];
				dst.score = src._score;
				dst.key = _entries[src._entryIndex]._key;
			}

			return resultsWithKey;
		}

	};

}

#endif // FUZZY_SEARCH_DATABASE_HPP
