// black 0,  white 255, a => 0 transparent, 255 opaque
#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CS248 {
vector<unsigned char> supersample_target; 

// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function
  // check bounds
  // cout << "Calling fill_sample" << endl;
	if (sx < 0 || sx >= target_w*this->sample_rate) return;
	if (sy < 0 || sy >= target_h*this->sample_rate) return;

  Color pixel_color;
	float inv255 = 1.0 / 255.0;
	pixel_color.r = supersample_target[4 * (sx + sy * target_w*this->sample_rate)] * inv255;
	pixel_color.g = supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 1] * inv255;
	pixel_color.b = supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 2] * inv255;
	pixel_color.a = supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 3] * inv255;

	pixel_color = ref->alpha_blending_helper(pixel_color, color);

	supersample_target[4 * (sx + sy * target_w*this->sample_rate)] = (uint8_t)(pixel_color.r * 255);
	supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 1] = (uint8_t)(pixel_color.g * 255);
	supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 2] = (uint8_t)(pixel_color.b * 255);
	supersample_target[4 * (sx + sy * target_w*this->sample_rate) + 3] = (uint8_t)(pixel_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= target_w) return;
	if (y < 0 || y >= target_h) return;

	// Color pixel_color;
	// float inv255 = 1.0 / 255.0;
	// pixel_color.r = render_target[4 * (x + y * target_w)] * inv255;
	// pixel_color.g = render_target[4 * (x + y * target_w) + 1] * inv255;
	// pixel_color.b = render_target[4 * (x + y * target_w) + 2] * inv255;
	// pixel_color.a = render_target[4 * (x + y * target_w) + 3] * inv255;

	// pixel_color = ref->alpha_blending_helper(pixel_color, color);

	// render_target[4 * (x + y * target_w)] = (uint8_t)(pixel_color.r * 255);
	// render_target[4 * (x + y * target_w) + 1] = (uint8_t)(pixel_color.g * 255);
	// render_target[4 * (x + y * target_w) + 2] = (uint8_t)(pixel_color.b * 255);
	// render_target[4 * (x + y * target_w) + 3] = (uint8_t)(pixel_color.a * 255);



  render_target[4 * (x + y * target_w)] = (uint8_t)(color.r * 255);
	render_target[4 * (x + y * target_w) + 1] = (uint8_t)(color.g * 255);
	render_target[4 * (x + y * target_w) + 2] = (uint8_t)(color.b * 255);
	render_target[4 * (x + y * target_w) + 3] = (uint8_t)(color.a * 255);


  for(int i=0; i<sample_rate; i++){
    for(int j=0; j<sample_rate; j++){
      fill_sample(x*sample_rate+i,y*sample_rate+j,color);
    }
  }
}

void SoftwareRendererImp::draw_svg( SVG& svg ) {
  printf("Inside draw_svg");
  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);

  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);  

  // resolve and send to render target
  resolve();


  unsigned char maxColor_r=0, maxColor_g=0, maxColor_b=0, maxColor_a=0, minColor_r=0, minColor_g=0, minColor_b=0, minColor_a=0;
  for(int x=0; x<this->target_w; x++){
    for(int y=0; y<this->target_h; y++){

      maxColor_r = max(maxColor_r, render_target[4 * (x + y * target_w)]);
      maxColor_g = max(maxColor_g, render_target[4 * (x + y * target_w) + 1]);
      maxColor_b = max(maxColor_b, render_target[4 * (x + y * target_w) + 2]);
      maxColor_a = max(maxColor_a, render_target[4 * (x + y * target_w) + 3]);

      minColor_r = min(minColor_r, render_target[4 * (x + y * target_w)]);
      minColor_g = min(minColor_g, render_target[4 * (x + y * target_w) + 1]);
      minColor_b = min(minColor_b, render_target[4 * (x + y * target_w) + 2]);
      minColor_a = min(minColor_a, render_target[4 * (x + y * target_w) + 3]);
    }
  }
  cout<<"\n";
  cout << "maxColor_r " << (int)maxColor_r << endl;
  cout << "maxColor_g " << (int)maxColor_g << endl;
  cout << "maxColor_b " << (int)maxColor_b << endl;
  cout << "maxColor_a " << (int)maxColor_a << endl;

  cout << "minColor_r " << (int)minColor_r << endl;
  cout << "minColor_g " << (int)minColor_g << endl;
  cout << "minColor_b " << (int)minColor_b << endl;
  cout << "minColor_a " << (int)minColor_a << endl;

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  cout<<"Setting sample rate to "<<sample_rate << endl;
  this->sample_rate = sample_rate;
  // supersample_target.resize(4*this->target_h*sample_rate*this->target_w*sample_rate);
  supersample_target.resize(4*this->target_h*4*this->target_w*4);
  std::fill(supersample_target.begin(),supersample_target.end(),255.0);
  cout<<"size of supersample_target (" << this->target_h*sample_rate<<", "<<this->target_w*sample_rate<<" ) =>"<<
  supersample_target.size()<<endl;
  // printf("Inside set_sample_rate");
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  cout << "Calling set_render_target"<< endl;
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  supersample_target.resize(4*height*sample_rate*width*sample_rate);
  std::fill(supersample_target.begin(),supersample_target.end(),255.0);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack
  transformation = transformation* element->transform;
	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}
  transformation = transformation* element->transform.inv();
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= target_w) return;
  if (sy < 0 || sy >= target_h) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  // render_target[4 * (sx + sy * target_w)] = (uint8_t)(color.r * 255);
  // render_target[4 * (sx + sy * target_w) + 1] = (uint8_t)(color.g * 255);
  // render_target[4 * (sx + sy * target_w) + 2] = (uint8_t)(color.b * 255);
  // render_target[4 * (sx + sy * target_w) + 3] = (uint8_t)(color.a * 255);

  // for(int i=0; i<sample_rate; i++){
  //   for(int j=0; j<sample_rate; j++){
  //     fill_sample(x*sample_rate+i,y*sample_rate+j,color);
  //   }
  // }
  fill_pixel(x, y, color);
}
void SoftwareRendererImp::line_helper_1(bool flipped,
                  float x0, float y0,
                  float x1, float y1,
                  Color color){
  int dx  = x1 - x0,
  dy  = y1 - y0,
  y   = y0,
  eps = 0;
  for ( int x = x0; x <= x1; x++ )  {
    if (flipped) rasterize_point(y,x,color);
    else rasterize_point(x,y,color);
    eps += dy;
    if ( (eps << 1) >= dx )  {
      y++;  eps -= dx;
    }
  }
}

void SoftwareRendererImp::line_helper_2(bool flipped,
                  float x0, float y0,
                  float x1, float y1,
                  Color color){
  int dx  = x1 - x0,
  dy  = y1 - y0,
  y   = y0,
  eps = 0;
  for ( int x = x0; x <= x1; x++ )  {
    if (flipped) rasterize_point(y,x,color);
    else rasterize_point(x,y,color);
    eps += dy;
    if ( (eps << 1) <= -dx )  {
      y--;  eps += dx;
    }
  }
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  ref->rasterize_line_helper(x0, y0, x1, y1, target_w, target_h, color, this);

  // float m = (y1-y0)/(x1-x0);
  // if(m>=0){
  //   if(x1<x0){//order points st. pt0 is to the left of pt1
  //     // swap(x0,x1); swap(y0,y1);
  //     if(m<1) line_helper_1(false, x1,y1,x0,y0,color); //increment x
  //     else line_helper_1(true, y1,x1,y0,x0,color); //increment y  
  //   }
  //   if(m<1) line_helper_1(false, x0,y0,x1,y1,color); //increment x
  //   else line_helper_1(true, y0,x0,y1,x1,color); //increment y
  // }
  // else{

  //   if(m>-1){
  //     if(x1<x0){//order points st. pt0 is to the left of pt1
  //       // swap(x0,x1); swap(y0,y1);
  //       line_helper_2(false, x1,y1,x0,y0,color); //increment x
  //     }
  //     line_helper_2(false, x0,y0,x1,y1,color); //increment x
  //   }
  //   else {
  //     if(y1<y0){//order points st. pt0 is below pt1
  //       // swap(x0,x1); swap(y0,y1);
  //       line_helper_2(true, y1,x1,y0,x0,color); //increment y
  //     }
  //     line_helper_2(true, y0,x0,y1,x1,color); //increment y

  //   }
  // }
  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

bool pointInTriangle(float* A, float* B, float*C, bool windDir, float x, float y){
  if(windDir){
    for(int i=0; i<3; i++){
      if((A[i]*x - B[i]*y + C[i]) < 0){
        return false;
      }
    }
  }
  else{
    for(int i=0; i<3; i++){
      if((A[i]*x - B[i]*y + C[i]) > 0){
        return false;
      }
    }
  }
  return true;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // return;
  // Task 1: 
  // Implement triangle rasterization
  float A[3], B[3], C[3];
  A[0] = y1-y0;
  A[1] = y2-y1;
  A[2] = y0-y2;

  B[0] = x1-x0;
  B[1] = x2-x1;
  B[2] = x0-x2;

  C[0] = y0*B[0] - x0*A[0];
  C[1] = y1*B[1] - x1*A[1];
  C[2] = y2*B[2] - x2*A[2];

  bool windDir = (A[0]*x2 - B[0]*y2 + C[0]) > 0;

  int minX = floor(min(x0, min(x1, x2)));
  int maxX = floor(max(x0, max(x1, x2)))+1;
  int minY = floor(min(y0, min(y1, y2)));
  int maxY = floor(max(y0, max(y1, y2)))+1;

  for (int x = minX; x<=maxX; x++){
    for (int y = minY; y<=maxY; y++){
      for(int i=0; i<sample_rate; i++){
        for(int j=0; j<sample_rate; j++){
      
          float cx = x+(0.5/sample_rate)+i/sample_rate, cy = y+(0.5/sample_rate)+j/sample_rate;
          if(pointInTriangle(A, B, C, windDir, cx,cy)) {
            fill_sample(x*sample_rate+i,y*sample_rate+j,color);
          }
        }
      }
    }
  }

  // Advanced Task
  // Implementing Triangle Edge Rules

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {
  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
  Color maxColor, minColor;
  for(int x=0; x<this->target_w; x++){
    for(int y=0; y<this->target_h; y++){
      // cout << "Calling resolve"<< endl;
      Color avg_color(0.0,0.0, 0.0, 0.0);
      for(int i=0; i<sample_rate; i++){
        for(int j=0; j< sample_rate; j++){
          int sx = x*this->sample_rate+i,  sy = y*this->sample_rate+j;
          avg_color.r += supersample_target[4*(sx+sy*target_w*sample_rate)]/255.0;// * supersample_target[4*(sx+sy*target_w*sample_rate)+3];
          avg_color.g += supersample_target[4*(sx+sy*target_w*sample_rate)+1]/255.0;// * supersample_target[4*(sx+sy*target_w*sample_rate)+3]; 
          avg_color.b += supersample_target[4*(sx+sy*target_w*sample_rate)+2]/255.0;// * supersample_target[4*(sx+sy*target_w*sample_rate)+3];
          avg_color.a += supersample_target[4*(sx+sy*target_w*sample_rate)+3]/255.0;  
        }
      }
      avg_color.r /= sample_rate*sample_rate;
      avg_color.g /= sample_rate*sample_rate;
      avg_color.b /= sample_rate*sample_rate;
      avg_color.a /= sample_rate*sample_rate;

      maxColor.r = max(maxColor.r, avg_color.r);
      maxColor.g = max(maxColor.g, avg_color.g);
      maxColor.b = max(maxColor.b, avg_color.b);
      maxColor.a = max(maxColor.a, avg_color.a);

      minColor.r = min(minColor.r, avg_color.r);
      minColor.g = min(minColor.g, avg_color.g);
      minColor.b = min(minColor.b, avg_color.b);
      minColor.a = min(minColor.a, avg_color.a);

      // render_target[4 * (x + y * target_w)] = (uint8_t)(avg_color.r * 255);
      // render_target[4 * (x + y * target_w) + 1] = (uint8_t)(avg_color.g * 255);
	    // render_target[4 * (x + y * target_w) + 2] = (uint8_t)(avg_color.b * 255);
	    // render_target[4 * (x + y * target_w) + 3] = (uint8_t)(avg_color.a * 255);

      fill_pixel(x,y,avg_color);
    }
  }
  cout<<"\n";
  cout << "maxColor.r " << maxColor.r << endl;
  cout << "maxColor.g " << maxColor.g << endl;
  cout << "maxColor.b " << maxColor.b << endl;
  cout << "maxColor.a " << maxColor.a << endl;

  cout << "minColor.r " << minColor.r << endl;
  cout << "minColor.g " << minColor.g << endl;
  cout << "minColor.b " << minColor.b << endl;
  cout << "minColor.a " << minColor.a << endl;
  
  // supersample_target.clear(); //free buffer of supersamples
  return;

}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248
