#include <iostream>
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

inline int coord_to_ind(const double, const double, const double, 
                        const int);
inline bool is_frame(const int, const int, const int, const int, 
                     const int, const int, const int);
inline void print_time(const chrono::system_clock::time_point,
                       const chrono::system_clock::time_point);


int main()
{
const int n_grid = 1000;  // pixels per axis
const int max_iter = 100; // maximum # of iterations
double x, y;              // coordinates of the complex plane
const double xmin = -2, xmax = 0.5;
const double ymin = -1.3, ymax = 1.3;
string filename = "./out/data.txt"; // name of the output file

// Define variables for printing a zoom window frame:
const int frame_width = 1; // frame width in pixels
const double xmin_frame = -0.4, xmax_frame = 0.2;
const double ymin_frame = 0.55, ymax_frame = 1.15;
const int m1 = coord_to_ind(xmin_frame, xmin, xmax, n_grid);
const int m2 = coord_to_ind(xmax_frame, xmin, xmax, n_grid);
const int n1 = coord_to_ind(ymin_frame, ymin, ymax, n_grid);
const int n2 = coord_to_ind(ymax_frame, ymin, ymax, n_grid);

// Open the output file:
ofstream outfile;
outfile.open(filename);
outfile << n_grid << " " << max_iter << " " << endl;
outfile << xmin << " " << xmax << " " << ymin << " " << ymax 
        << " " << endl;

// Set the start time of the Mandelbrot set algorithm:
chrono::system_clock::time_point t1 = chrono::system_clock::now();


// Calculate the Mandelbrot set.
// Go through the complex plane and iterate the function 
// f(z,c) = z^2 + c at z = 0 for each value of c = x + i*y.
int cnt_iter; // # of iterations
double zr, zi, zr_old;

for (int m = 0; m != n_grid; ++m)
{
  x = xmin + (xmax-xmin)*m/(n_grid-1);

  for (int n = 0; n != n_grid; ++n)
  {
    y = ymin + (ymax-ymin)*n/(n_grid-1);

    // f(z,c) = z^2 + c
    // c = x + i*y,  z = zr + i*zi
    for (cnt_iter = 0, zr = 0, zi = 0; 
         zr*zr + zi*zi < 4 && cnt_iter != max_iter; 
         ++cnt_iter)
    {
      zr_old = zr;
      zr = zr*zr - zi*zi + x;
      zi = 2*zr_old*zi + y;
    }

    // Set the value of the frame pixels to nan:
    if(is_frame(m, n, m1, m2, n1, n2, frame_width))
      outfile << "nan" << " ";
    // Set the pixels in the Mandelbrot set to zero:
    else if (cnt_iter == max_iter)
      outfile << 0 << " ";
    // and the other pixels to the # of iterations:
    else
      outfile << cnt_iter << " ";
  }
}
outfile << endl;
outfile.close();

// Set the end time:
chrono::system_clock::time_point t2 = chrono::system_clock::now();

cout << filename << endl;
cout << "n_grid = " << n_grid << ", max_iter = " << max_iter << endl;
print_time(t1, t2);

return 0;
}



//  For a given coordinate value, return the corresponding index.
inline int coord_to_ind(const double val, const double minv, const double maxv, 
                        const int ngrid)
{
  return (val-minv)*(ngrid-1)/(maxv-minv);
}


// Check if a pixel is part of the zoom window frame.
inline bool is_frame(const int m, const int n, const int m1, const int m2, 
                     const int n1, const int n2, const int w)
{
  return 
  (
  ((m >= m1 && m <= m1+w) || (m >= m2 && m <= m2+w)) && (n >= n1 && n <= n2)
  ) || (
  (m >= m1 && m <= m2) && ((n >= n1 && n <= n1+w) || (n >= n2 && n <= n2+w))
  );
}


// Print the execution time.
inline void print_time(const chrono::system_clock::time_point t1,
                       const chrono::system_clock::time_point t2)
{
  auto duration = chrono::duration_cast<chrono::seconds>(t2 - t1).count();
  int hours = duration/3600;
  int minutes = duration/60 - hours*60;
  int seconds = duration - hours*3600 - minutes*60;
  cout << hours << " h " << minutes << " min " << seconds << " s" << endl;
}
