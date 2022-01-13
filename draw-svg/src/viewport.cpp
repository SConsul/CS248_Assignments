#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvans coordinates.
  this->x = x;
  this->y = y;
  this->span = span; 

  double data[9] = {0.5/span, 0.0, 0.5-0.5*x/span, 0.0, 0.5/span, 0.5-0.5*y/span, 0.0, 0.0, 1.0};
  Matrix3x3 m(data);
  set_canvas_to_norm(m);

  
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
