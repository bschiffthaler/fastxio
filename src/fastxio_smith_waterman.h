#ifndef _FASTX_IO_SW_H_
#define _FASTX_IO_SW_H_

#include <string>
#include <fstream>
#include <vector>
#include <matrix.h>

namespace FASTX {

  /**
   * @brief Calculate the highest scoring neighbor in the SW matrix
   *
   * To find the optimal path in reverse from the maximal 
   * overall value of the matrix. In case of ties, the 
   * function prioritizes the upper left diagonal, then the
   * upper row, and lastly the left column.
   *
   * @param a The upper left diagonal value
   * @param b The upper row value
   * @param c The left column value
   * 
   * @return The value of the highest scorng neighbor
   */
  score_t sw_max_neighbour(score_t a, score_t b, score_t c);

  /**
   * @brief Helper struct to store the matrix index in the Smith
   * Waterman scoring matrix
   */
  struct SW_path
  {
    unsigned long i = 0; /**< Row index in matrix */
    unsigned long j = 0; /**< Column index in matrix */
  };

  /** 
   * @brief A Smith-Waterman class to perform local alignments.
   *
   * This class should be treated as a very simple, early
   * version of an SW aligner. It performs very basic local
   * alignments of short strings. For long strings, memory
   * efficiency should be a concern as it builds a scoring
   * matrix (Na * Nb), where Na is the length of string A,
   * and Nb is the length of string B.
   *
   */
  class SW {
  public:
    /**
     * @brief Perform a Smith-Waterman alignment of two strings
     *
     * @param s1 A DNA string
     * @param s2 Another DNA string
     */
    SW(std::string s1, std::string s2);

    /**
     * @brief Print the scoring matrix
     *
     * @param out The sink, derived from `std::ostream`
     */
    void print_matrix(std::ostream& out);

    /**
     * @brief Print the path that was calculated through the matrix
     *
     * @param out The sink, derived from `std::ostream`
     */
    void print_path(std::ostream& out);

    /**
     * @brief Print the local alignment
     *
     * @param out The sink, derived from `std::ostream`
     * @param sw The Smith-Waterman class
     */
    friend std::ostream& operator<<(std::ostream& out, const SW sw);

  private:
    std::string _s1;
    std::string _s2;
    Matrix<score_t> _scores;
    std::vector<SW_path> _path;
    unsigned long _s1_start;
    unsigned long _s1_stop;
    unsigned long _s2_start;
    unsigned long _s2_stop;
  };

  /**
   * @brief Print the local alignment
   *
   * @param out The sink, derived from `std::ostream`
   * @param sw The Smith-Waterman class
   */
  std::ostream& operator<<(std::ostream& out, const SW sw);
  
}
#endif
