# ifndef HEADER_H_
# define HEADER_H_

# include <iostream>
# include <vector>
# include <algorithm>
# include <tuple>
# include <utility>
# include <stdexcept>
# include <numeric>
# include <type_traits>
# include <concepts>
# include <random>
# include <variant>
# include <cmath>
# include <fstream>
# include <memory>
# include <array>
# include <limits>
# include <string>
# include <sstream>
# include <chrono>
# include "timer.h"

# define USE_TIMER true

# if USE_TIMER
# define TIME_THIS() Timer funcTimer(__func__)
# else
# define TIME_THIS()
# endif


# endif 