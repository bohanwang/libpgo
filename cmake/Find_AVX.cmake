# https://gist.github.com/UnaNancyOwen/263c243ae1e05a2f9d0e
# Check for the presence of AVX and figure out the flags to use for it.
macro(CHECK_FOR_AVX)
  include(CheckCXXSourceRuns)

  # Check AVX
  if(MSVC AND NOT MSVC_VERSION LESS 1600)
    set(CMAKE_REQUIRED_FLAGS "/arch:AVX")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    set(CMAKE_REQUIRED_FLAGS "-march=native -mtune=native")
  else()
    message(FATAL_ERROR "Unsupported compiler")
  endif()

  check_cxx_source_runs("
    #include <immintrin.h>
    int main()
    {
      __m256 a, b, c;
      const float src[8] = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
      float dst[8];
      a = _mm256_loadu_ps( src );
      b = _mm256_loadu_ps( src );
      c = _mm256_add_ps( a, b );
      _mm256_storeu_ps( dst, c );

      for( int i = 0; i < 8; i++ ){
        if( ( src[i] + src[i] ) != dst[i] ){
          return -1;
        }
      }

      return 0;
    }"
  HAVE_AVX_EXTENSIONS)

  # Check AVX2
  if(MSVC AND NOT MSVC_VERSION LESS 1800)
    set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    set(CMAKE_REQUIRED_FLAGS "-march=native -mtune=native")
  else()
    message(FATAL_ERROR "Unsupported compiler")
  endif()

  check_cxx_source_runs("
    #include <immintrin.h>
    int main()
    {
      __m256i a, b, c;
      const int src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
      int dst[8];
      a =  _mm256_loadu_si256( (__m256i*)src );
      b =  _mm256_loadu_si256( (__m256i*)src );
      c = _mm256_add_epi32( a, b );
      _mm256_storeu_si256( (__m256i*)dst, c );

      for( int i = 0; i < 8; i++ ){
        if( ( src[i] + src[i] ) != dst[i] ){
          return -1;
        }
      }

      return 0;
    }"
  HAVE_AVX2_EXTENSIONS)

  # from pytorch
  check_cxx_source_runs("
    #include <immintrin.h>
  
    int main()
    {
      __m512i a = _mm512_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0);
      __m512i b = a;
      __mmask64 equality_mask = _mm512_cmp_epi8_mask(a, b, _MM_CMPINT_EQ);
      return 0;
    }"
  HAVE_AVX512_EXTENSIONS)

endmacro(CHECK_FOR_AVX)
