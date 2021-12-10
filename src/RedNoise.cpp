#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

#define WIDTH 640
#define HEIGHT 480

bool is_number(const std::string &s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

std::vector<CanvasPoint> interpolatePoints(CanvasPoint left, CanvasPoint right, int steps) {	
	float stepX = (right.x - left.x) / (steps - 1);
	float stepY = (right.y - left.y) / (steps - 1);
	float stepZ = (right.depth - left.depth) / (steps - 1);

	CanvasPoint curr = left;
	std::vector<CanvasPoint> v;
	v.push_back(curr);

	for (int i = 0; i < steps - 1; ++i) {
		curr.x = curr.x + stepX;
		curr.y = curr.y + stepY;
		curr.depth = curr.depth + stepZ;
		v.push_back(curr);
	}

	return v;
}

std::vector<TexturePoint> interpolatePoints(TexturePoint left, TexturePoint right, int steps) {
	float stepX = (right.x - left.x) / (steps - 1);
	float stepY = (right.y - left.y) / (steps - 1);

	TexturePoint curr = left;
	std::vector<TexturePoint> v;
	v.push_back(curr);

	for (int i = 0; i < steps - 1; ++i) {
		curr.x = curr.x + stepX;
		curr.y = curr.y + stepY;
		v.push_back(curr);
	}

	return v;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, std::vector <std::vector <float>> &depthBuffer) {
	float xDiff = to.x - from.x;
 	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
 	float numberOfSteps = fmax(fabs(xDiff), fmax(fabs(yDiff),fabs(zDiff)));
 	float xStepSize = xDiff/numberOfSteps;
 	float yStepSize = yDiff/numberOfSteps;
	float zStepSize = zDiff/numberOfSteps;

 	for (float i=0.0; i<numberOfSteps; i++) {
   		float x = from.x + (xStepSize*i);  //int
   		float y = from.y + (yStepSize*i); // int
		float z = -(from.depth + (zStepSize*i));
		uint32_t c = (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
		if (z == 0) window.setPixelColour(round(x), round(y), c);
			else if (depthBuffer[round(x)][round(y)] < 1.0/z) {
				window.setPixelColour(round(x), round(y), c);
				depthBuffer[round(x)][round(y)] = 1.0/z;
			}
	 }
}

void drawTriangle(DrawingWindow &window, CanvasTriangle tr, Colour colour, std::vector <std::vector <float>> &depthBuffer) {
	CanvasPoint x = tr.vertices[0];
	CanvasPoint y = tr.vertices[1];
	CanvasPoint z = tr.vertices[2];

	drawLine(window, x, y, colour, depthBuffer);
	drawLine(window, x, z, colour, depthBuffer);
	drawLine(window, y, z, colour, depthBuffer);
}

void fillTopTriangle(DrawingWindow &window, CanvasTriangle tr, Colour colour, std::vector <std::vector <float>> &depthBuffer) {
	float invslope1 = (tr.vertices[1].x - tr.vertices[0].x) / (tr.vertices[1].y - tr.vertices[0].y);
	float invslope2 = (tr.vertices[2].x - tr.vertices[0].x) / (tr.vertices[2].y - tr.vertices[0].y);
	float invslope3 = (tr.vertices[1].depth - tr.vertices[0].depth) / (tr.vertices[1].y - tr.vertices[0].y);
	float invslope4 = (tr.vertices[2].depth - tr.vertices[0].depth) / (tr.vertices[2].y - tr.vertices[0].y);
	float curx1 = tr.vertices[0].x;
	float curx2 = tr.vertices[0].x;
	float curz1 = tr.vertices[0].depth;
	float curz2 = tr.vertices[0].depth;
	for (float i = tr.vertices[0].y; i < tr.vertices[1].y; ++i) {
		drawLine(window, CanvasPoint(curx1, i, curz1), CanvasPoint(curx2, i, curz2), colour, depthBuffer);
		curx1 += invslope1;
		curx2 += invslope2;
		curz1 += invslope3;
		curz2 += invslope4;
	}
}

void fillBottomTriangle(DrawingWindow &window, CanvasTriangle tr, Colour colour, std::vector <std::vector <float>> &depthBuffer) {
	float invslope1 = (tr.vertices[2].x - tr.vertices[0].x) / (tr.vertices[2].y - tr.vertices[0].y);
	float invslope2 = (tr.vertices[2].x - tr.vertices[1].x) / (tr.vertices[2].y - tr.vertices[1].y);
	float invslope3 = (tr.vertices[2].depth - tr.vertices[0].depth) / (tr.vertices[2].y - tr.vertices[0].y);
	float invslope4 = (tr.vertices[2].depth - tr.vertices[1].depth) / (tr.vertices[2].y - tr.vertices[1].y);
  	float curx1 = tr.vertices[2].x;
  	float curx2 = tr.vertices[2].x;
	float curz1 = tr.vertices[2].depth;
	float curz2 = tr.vertices[2].depth;
	for (float i = tr.vertices[2].y; i > tr.vertices[1].y; --i) {
		drawLine(window, CanvasPoint(curx1, i, curz1), CanvasPoint(curx2, i, curz2), colour, depthBuffer);
		curx1 -= invslope1;
		curx2 -= invslope2;
		curz1 -= invslope3;
		curz2 -= invslope4;
	}
}

void rasterTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour c, std::vector <std::vector <float>> &depthBuffer) {
	
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (triangle.vertices[i].y > triangle.vertices[j].y) std::swap(triangle.vertices[i], triangle.vertices[j]);
			
	drawTriangle(window, triangle, c, depthBuffer);	

	float yDist = triangle.vertices[0].y - triangle.vertices[2].y;
	float xDist = triangle.vertices[0].x - triangle.vertices[2].x;
	float zDist = triangle.vertices[0].depth - triangle.vertices[2].depth;
	
	float ratio = xDist / yDist;
	float ratio2 = zDist / yDist;
	
	float yDist2 = triangle.vertices[1].y - triangle.vertices[2].y;

	float new_x = triangle.vertices[2].x + (yDist2 * ratio);
	float new_z = triangle.vertices[2].depth + (yDist2 * ratio2);
	
	CanvasPoint t = CanvasPoint(new_x, triangle.vertices[1].y, new_z);
	drawLine(window, t, triangle.vertices[1], c, depthBuffer);		

	CanvasTriangle topTr = CanvasTriangle(triangle.vertices[2], triangle.vertices[1], t);
	CanvasTriangle botTr = CanvasTriangle(triangle.vertices[1],t, triangle.vertices[0]);
	//drawTriangle(window, topTr, c);
	//drawTriangle(window, botTr, c);
	fillTopTriangle(window, topTr, c, depthBuffer);
	fillBottomTriangle(window, botTr, c, depthBuffer);	
}

void rasterTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture, std::vector <std::vector <float>> &depthBuffer) {

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (triangle.vertices[i].y > triangle.vertices[j].y) std::swap(triangle.vertices[i], triangle.vertices[j]);

	
	float yDist = triangle.vertices[0].y - triangle.vertices[2].y;
	float xDist = triangle.vertices[0].x - triangle.vertices[2].x;
	float zDist = triangle.vertices[0].depth - triangle.vertices[2].depth;

	
	float ratio = xDist / yDist;
	float ratio2 = zDist / yDist;

	float yDist2 = triangle.vertices[1].y - triangle.vertices[2].y;
	

	float new_x = triangle.vertices[2].x + (yDist2 * ratio);
	float new_z = triangle.vertices[2].depth + (yDist2 * ratio2);

	CanvasPoint t = CanvasPoint(new_x, triangle.vertices[1].y, new_z);

	float ratio_t = (triangle.vertices[1].y - triangle.vertices[2].y) / (triangle.vertices[0].y - triangle.vertices[2].y);
	t.texturePoint.x = triangle.vertices[2].texturePoint.x + ratio_t * (triangle.vertices[0].texturePoint.x - triangle.vertices[2].texturePoint.x);
	t.texturePoint.y = triangle.vertices[2].texturePoint.y + ratio_t *  (triangle.vertices[0].texturePoint.y - triangle.vertices[2].texturePoint.y);

	//TOP triangle
	int step = triangle.vertices[1].y - triangle.vertices[2].y + 2;
	std::vector<CanvasPoint> l = interpolatePoints(triangle.vertices[2], triangle.vertices[1], step);
	std::vector<CanvasPoint> r = interpolatePoints(triangle.vertices[2], t, step);
	std::vector<TexturePoint> l_t = interpolatePoints(triangle.vertices[2].texturePoint, triangle.vertices[1].texturePoint, step);
	std::vector<TexturePoint> r_t = interpolatePoints(triangle.vertices[2].texturePoint, t.texturePoint, step);

	for (int i = 0; i < l.size(); ++i) {
		int step = abs(l[i].x - r[i].x);
		
		std::vector<CanvasPoint> p = interpolatePoints(l[i], r[i], step + 2);
		std::vector<TexturePoint> p_t = interpolatePoints(l_t[i], r_t[i], step + 2);
		
		for (int j = 0; j < p_t.size(); ++j) {
			int x = round(p[j].x);
			int y = round(p[i].y);
			float z = -p[j].depth;
			uint32_t c = texture.pixels[int(p_t[j].x) + int(p_t[j].y) * texture.width];
			if (z == 0) window.setPixelColour(x, y, c);
				else if (depthBuffer[x][y] < 1.0/z) {
					window.setPixelColour(x, y, c);
					depthBuffer[x][y] = 1.0/z;
				}
			//window.setPixelColour(x, y, c);
		}

	}
	
	step = triangle.vertices[0].y - triangle.vertices[1].y + 2;
	l = interpolatePoints(triangle.vertices[0], triangle.vertices[1], step);
	r = interpolatePoints(triangle.vertices[0], t, step);
	l_t = interpolatePoints(triangle.vertices[0].texturePoint, triangle.vertices[1].texturePoint, step);
	r_t = interpolatePoints(triangle.vertices[0].texturePoint, t.texturePoint, step);

	for (int i = 0; i < l.size(); ++i) {
		step = abs(l[i].x - r[i].x);
		
		std::vector<CanvasPoint> p = interpolatePoints(l[i], r[i], step + 2);
		std::vector<TexturePoint> p_t = interpolatePoints(l_t[i], r_t[i], step + 2);
		
		for (int j = 0; j < p_t.size(); ++j) {
			int x = round(p[j].x);
			int y = round(p[i].y);
			float z = -p[j].depth;
			uint32_t c = texture.pixels[int(p_t[j].x) + int(p_t[j].y) * texture.width];
			if (z == 0) window.setPixelColour(x, y, c);
				else if (depthBuffer[x][y] < 1.0/z) {
					window.setPixelColour(x, y, c);
					depthBuffer[x][y] = 1.0/z;
				}
			//window.setPixelColour(x, y, c);
		}
		
	}
	
}

std::unordered_map <std::string, Colour> parseMtl(std::string filename, std::unordered_map <std::string, TextureMap> &textures) {
	std::unordered_map <std::string, Colour> colours;

	std::ifstream in(filename);
	if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        //exit(1);
    }

	std::string line;
	std::string name;
	while(std::getline(in, line)) {
		std::vector <std::string> x = split(line, ' ');
		if (x[0] == "newmtl") {
			name = x[1];
		}
		if (x[0] == "Kd") {
			int R = (int) (stof(x[1]) * 255.0);
			int G = (int) (stof(x[2]) * 255.0);
			int B = (int) (stof(x[3]) * 255.0);
			Colour colour = Colour(name, R, G, B);
			colours[name] = colour;
		}

		if (x[0] == "map_Kd") {
			textures[name] = TextureMap(x[1]);
		}
	}

	return colours;
}

std::vector<ModelTriangle> parseObj(std::string filename, float scale_factor, std::unordered_map <std::string, TextureMap> &textures, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, TextureMap MAP) {

	std::ifstream in(filename);
  	if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        exit(1);
    }

	vertexNormals.push_back(glm::vec3(0,0,0));

	std::string line;
	std::vector <ModelTriangle> arr;
	std::vector <glm::vec3> v;
	std::vector <TexturePoint> texturePoints;
	texturePoints.push_back(TexturePoint(0, 0));
	v.push_back(glm::vec3(0, 0, 0));
	std::unordered_map <std::string, Colour> colours; 
	std::string colour;
	while(std::getline(in, line)) {
		std::vector <std::string> x = split(line, ' ');
		if (x[0] == "mtllib") {
			colours = parseMtl(x[1], textures);
		}
		if (x[0] == "usemtl") colour = x[1];
		if (x[0] == "v") {
			glm::vec3 vertex(stof(x[1]) * scale_factor, stof(x[2]) * scale_factor, stof(x[3]) * scale_factor);
			vertexNormals.push_back(glm::vec3(0,0,0));
			v.push_back(vertex);
		}
		if (x[0] == "vt") {
			TexturePoint vertex(stof(x[1])* MAP.width, stof(x[2]) * MAP.height);
			texturePoints.push_back(vertex);
		}
		if (x[0] == "f") {
			int ind[4] = {0, 0, 0, 0};
			int t_ind[4] = {0, 0, 0, 0};

			for (int i = 1; i <= 3; ++i) {
				std::vector <std::string> vals = split(x[i], '/');
				if (is_number(vals[1])) t_ind[i] = stoi(vals[1]);
				ind[i] = stoi(vals[0]);
			}
			ModelTriangle triangle = ModelTriangle(v[ind[1]], v[ind[2]], v[ind[3]], colours[colour]);
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			triangle.normal = glm::normalize(glm::cross(e1, e0)) * (-1.f);
			
			vertexNormals[ind[1]] += triangle.normal;
			vertexNormals[ind[2]] += triangle.normal;
			vertexNormals[ind[3]] += triangle.normal;
			
			std::array<TexturePoint, 3> texture_arr;
			for (int i = 1; i <= 3; ++i) {
				texture_arr[i - 1] = texturePoints[t_ind[i]];
			}			
			triangle.texturePoints = texture_arr;
			
			std::vector <int> aux;
			aux.push_back(ind[1]);
			aux.push_back(ind[2]);
			aux.push_back(ind[3]);
			trToNormal.push_back(aux);
			arr.push_back(triangle);
		}
	}

	for (int i = 1; i < vertexNormals.size(); ++i)
		vertexNormals[i] = glm::normalize(vertexNormals[i]);


	return arr;			
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 &cameraPosition, glm::vec3 vertexPosition, float focalLength, glm::mat3 &cameraOrientation) {
	vertexPosition -= cameraPosition;
	vertexPosition = vertexPosition * cameraOrientation;
	float x = focalLength * vertexPosition[0] / vertexPosition[2] * (-500) + WIDTH / 2;
	float y = focalLength * vertexPosition[1] / vertexPosition[2] * 500 + HEIGHT / 2;
		
	return CanvasPoint(x, y, vertexPosition.z);
}

glm::vec3 get3DIntersectionPoint(CanvasPoint point, glm::vec3 &cameraPosition, float focalLength, glm::mat3 &cameraOrientation) {
	float x = point.x;
	float y = point.y;
	float z = point.depth;

	glm::vec3 vertexPosition((x - WIDTH / 2) / 1000, (y - HEIGHT / 2) / (-1000), z);
	vertexPosition = vertexPosition * glm::inverse(cameraOrientation);
	//vertexPosition -= cameraPosition;

	return vertexPosition;
}

/* week4
void task6(DrawingWindow &window, std::vector <ModelTriangle> &v, glm::vec3 camPos, float focalLength) {
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	for (auto &it: v) {
		for (int i = 0; i < 3; ++i)  {
			CanvasPoint point = getCanvasIntersectionPoint(camPos, it.vertices[i], focalLength);
			std::cout << point.x << " "<< point.y<<std::endl;
			window.setPixelColour(round(point.x ), round(point.y), colour);
		}
	}
}*/

void orbit(glm::vec3 &cameraPosition, float speed, bool zoom, bool sphere) {
	float rad = speed;
	glm::mat3 yRotation = glm::mat3(
		cos(rad), 0, -sin(rad),
		0, 1, 0,
		sin(rad), 0, cos(rad)
		);
	if (!sphere) {
	glm::vec3 distance = (cameraPosition - glm::vec3(0,0,0)) * 0.05f;
	if (!zoom) cameraPosition = cameraPosition - speed*distance;
		else cameraPosition = cameraPosition + speed*distance;
	cameraPosition = cameraPosition * yRotation;
	} else {
		cameraPosition = cameraPosition * yRotation;
	}
}

void lookAt(glm::mat3 &cameraOrientation, glm::vec3 point, glm::vec3 &cameraPosition) {
	glm::vec3 forward = glm::normalize(cameraPosition - point);
	glm::vec3 right = glm::cross(glm::normalize(glm::vec3(0, 1, 0)), forward);
	glm::vec3 up = glm::cross(forward, right);
	cameraOrientation = glm::mat3(right, up, forward);
 }

void drawWireframeScene(DrawingWindow &window, std::vector <ModelTriangle> &v, glm::vec3 &camPos, float focalLength, glm::mat3 &cameraOrientation, std::vector <std::vector <float>> &depthBuffer, float &rad, bool zoom) {
	window.clearPixels();
	orbit(camPos, rad, zoom, false);
	lookAt(cameraOrientation, glm::vec3(0,0,0), camPos);
	for (auto &it: depthBuffer) std::fill(it.begin(), it.end(), 0);

	for (auto &it: v) {
		CanvasPoint vertices[3];
		for (int i = 0; i < 3; ++i)  {
			vertices[i] = getCanvasIntersectionPoint(camPos, it.vertices[i], focalLength, cameraOrientation);
		}
		drawTriangle(window, CanvasTriangle(vertices[0], vertices[1], vertices[2]), Colour(255,255,255), depthBuffer);
	}
}

void drawRasterisedScene(DrawingWindow &window, std::vector <ModelTriangle> &v, glm::vec3 &cameraPosition, float focalLength, std::vector <std::vector <float>> &depthBuffer, glm::mat3 &cameraOrientation, std::unordered_map <std::string, TextureMap> &textures, float &rad, bool zoom) {
	window.clearPixels();
	orbit(cameraPosition, rad, zoom, false);
	lookAt(cameraOrientation, glm::vec3(0,0,0), cameraPosition);
	for (auto &it: depthBuffer) std::fill(it.begin(), it.end(), 0);
	
	for (auto &it: v) {
		CanvasPoint vertices[3];
		for (int i = 0; i < 3; ++i)  {
			vertices[i] = getCanvasIntersectionPoint(cameraPosition, it.vertices[i], focalLength, cameraOrientation);
		}
		if (it.texturePoints[0].x == 0 && it.texturePoints[0].y == 0 && it.texturePoints[1].x == 0 && it.texturePoints[1].y == 0 && it.texturePoints[2].x == 0 && it.texturePoints[2].y== 0) {
			rasterTriangle(window, CanvasTriangle(vertices[0], vertices[1], vertices[2]), it.colour, depthBuffer);
		}
		else {
			for (int i = 0; i < 3; ++i) {
				vertices[i].texturePoint.x = it.texturePoints[i].x;
				vertices[i].texturePoint.y = it.texturePoints[i].y;
			}

			rasterTexturedTriangle(window, CanvasTriangle(vertices[0], vertices[1], vertices[2]),  TextureMap(textures[it.colour.name]), depthBuffer);
		}
	}
}

RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 &rayDirection, std::vector<ModelTriangle> &v, int &noIntTriangles, glm::vec3 &barycentric) {
	RayTriangleIntersection best;
	best.distanceFromCamera = FLT_MAX;
	best.triangleIndex = -1;
	int triangleIndex = 0;

	for (auto &triangle: v) {
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution[0], u = possibleSolution[1], v = possibleSolution[2];

		if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) &&  ((u + v) <= 1.0) && (t >= 0.0)) {
			noIntTriangles++; 
			if (possibleSolution[0] < best.distanceFromCamera) {
				glm::vec3 intersectionPoint = cameraPosition + t * rayDirection;
				//glm::vec3 intersectionPoint = triangle.vertices[0] + u * e0 + v * e1;
				RayTriangleIntersection x = RayTriangleIntersection(intersectionPoint, t, triangle, triangleIndex);
				best = x;
				barycentric = glm::vec3(u, v, 1.0f - u - v);
			}
		}
		triangleIndex++;
	}

	return best;
}

uint32_t getColour(float red, float green, float blue, float diffuse, float specular, float ambientStrength) {
	red = fmin(red * (diffuse + specular + ambientStrength), 255.0);
	green = fmin(green * (diffuse + specular + ambientStrength), 255.0);
	blue = fmin(blue * (diffuse + specular + ambientStrength), 255.0);
	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	return colour;
}

float getDiffuse(float r, glm::vec3 normal, glm::vec3 lightDirection) {
	float angleDot = fmax(glm::dot(normal, lightDirection), 0.0); //angle of incidence lighting
	float brightness = 1.0 / (4.0 * 3.14 * r * r); // proximity lighting
	float diffuse = angleDot * brightness; // diffuse lighting
	return diffuse;
}

float getSpecular(glm::vec3 lightDirection, glm::vec3 normal, glm::vec3 &cameraPosition, glm::vec3 intersectionPoint, int specularExponent) {
	glm::vec3 Ri = (-1.0f) * lightDirection;
	glm::vec3 reflection = Ri - 2.0f * normal * (glm::dot(Ri, normal));
	glm::vec3 view = glm::normalize(cameraPosition - intersectionPoint);
	float specular = fmax(fmin(glm::dot(view, reflection), 1), 0.0f); // [0, 1]
	specular = pow(specular, specularExponent); // specular lighting
	return specular;
}

std::pair< glm::vec3, glm::vec2> getNormal(glm::vec3 barycentric, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, ModelTriangle triangle, std::vector<ModelTriangle> &v) {
	for (int k = 0; k < v.size(); ++k) {
		bool OK = true;
		for (int i = 0; i < 3; ++i)
			if (triangle.vertices[i] != v[k].vertices[i]) {
				OK = false;
				break;
			}
		//if i found the triangle
		if (OK) {
			std::vector<int> index = trToNormal[k];
			glm::vec3 normal1 = vertexNormals[index[0]] * barycentric[2];
			glm::vec3 normal2 = vertexNormals[index[1]] * barycentric[0];
			glm::vec3 normal3 = vertexNormals[index[2]] * barycentric[1];

			glm::vec2 t_p1 = glm::vec2(triangle.texturePoints[0].x, triangle.texturePoints[0].y) * barycentric[2];
			glm::vec2 t_p2 = glm::vec2(triangle.texturePoints[1].x, triangle.texturePoints[1].y) * barycentric[0];
			glm::vec2 t_p3 = glm::vec2(triangle.texturePoints[2].x, triangle.texturePoints[2].y) * barycentric[1];
		
			return {glm::normalize(normal1 + normal2 + normal3), (t_p1 + t_p2 + t_p3)};
		}
	}
	return {glm::vec3(0,0,0), glm::vec2(0,0)};
}

std::pair<Colour, std::pair<std::array<TexturePoint, 3>, glm::vec2> > mirror(glm::vec3 intersectionPoint, glm::vec3 rayDirection, std::vector<ModelTriangle> &v, int depth, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, glm::vec3 prevNormal, ModelTriangle &reflected_triangle) {
	glm::vec3 barycentric(0,0,0);
	int sample2 = 0;
	std::array<TexturePoint, 3> x;
	if (depth == 10) return {Colour(0,0,0), {x, glm::vec2(0,0)}};
	depth++;

	RayTriangleIntersection triangleInter = getClosestIntersection(intersectionPoint + prevNormal * 0.0001f, rayDirection, v, sample2, barycentric);
	if (triangleInter.triangleIndex == -1) return {Colour(0,0,0), {x, glm::vec2(0,0)}};
	reflected_triangle = triangleInter.intersectedTriangle;
	std::pair<glm::vec3, glm::vec2> NORMALES = getNormal(barycentric, vertexNormals, trToNormal, triangleInter.intersectedTriangle, v);
	glm::vec3 NORMAL = NORMALES.first;

	if (triangleInter.intersectedTriangle.colour.name == "Blue") {
		glm::vec3 Ri = rayDirection;
		glm::vec3 reflection = Ri - 2.0f * NORMAL * (glm::dot(Ri, NORMAL));
		
		return mirror(triangleInter.intersectionPoint, reflection, v, depth, vertexNormals, trToNormal, NORMAL, reflected_triangle);
	}
	return {triangleInter.intersectedTriangle.colour, {triangleInter.intersectedTriangle.texturePoints, NORMALES.second}};
}


glm::vec3 refract(glm::vec3 rayDirection, glm::vec3 NORMAL, float ior) {
  float h = glm::dot(-rayDirection, NORMAL), cos = 1.0f;
  if (h < 1) cos = h;
  glm::vec3 P = (rayDirection + cos * NORMAL) * ior;
  float r = float((-1) * std::sqrt(std::fabs(1.0 - glm::dot(P, P))));
  glm::vec3 PLL = r * NORMAL;
  return (P + PLL);
}

std::pair<Colour, std::pair<std::array<TexturePoint, 3>, glm::vec2> > refra(glm::vec3 intersectionPoint, glm::vec3 rayDirection, std::vector<ModelTriangle> &v, int depth, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, glm::vec3 prevNormal, ModelTriangle &reflected_triangle, bool inside) {
	glm::vec3 barycentric(0,0,0);
	int sample2 = 0;
	std::array<TexturePoint, 3> x;
	if (depth == 10) return {Colour(0,0,0), {x, glm::vec2(0,0)}};
	depth++;

	RayTriangleIntersection triangleInter = getClosestIntersection(intersectionPoint + prevNormal * 0.1f, rayDirection, v, sample2, barycentric);
	if (triangleInter.triangleIndex == -1) return {Colour(0,0,0), {x, glm::vec2(0,0)}};
	reflected_triangle = triangleInter.intersectedTriangle;
	std::pair<glm::vec3, glm::vec2> NORMALES = getNormal(barycentric, vertexNormals, trToNormal, triangleInter.intersectedTriangle, v);
	glm::vec3 NORMAL = NORMALES.first;
	if (triangleInter.intersectedTriangle.colour.name == "Red") {
		glm::vec3 refraction(0,0,0);
		if (!inside) {
          refraction = refract(rayDirection, NORMAL, 1.5f);
          inside = true;
      } else {
          NORMAL = NORMAL * -1.f;
          refraction = refract(rayDirection, NORMAL, 1/1.5f);
          inside = false;
      }
		return refra(triangleInter.intersectionPoint, refraction, v, depth, vertexNormals, trToNormal, NORMAL, reflected_triangle, inside);
	}
	return {triangleInter.intersectedTriangle.colour, {triangleInter.intersectedTriangle.texturePoints, NORMALES.second}};
}



void drawRayTracedScene(DrawingWindow &window, glm::vec3 &cameraPosition, float focalLength, glm::mat3 &cameraOrientation, std::vector<ModelTriangle> &v, glm::vec3 &lightLocation, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, TextureMap &MAP, std::vector <glm::vec3> &lights, bool REFLEC,  bool REFLEC2, bool zoom, bool SPHERE, bool rotate) {
	window.clearPixels();
	if (rotate) {
		orbit(cameraPosition, 0.05, zoom, SPHERE);
		lookAt(cameraOrientation, glm::vec3(0,0,0), cameraPosition);
	}
	int noIntTriangles = 0;
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			glm::vec3 barycentric(0, 0, 0);
			glm::vec3 empty(0, 0, 0);
 			
			glm::vec3 rayDirection = glm::normalize(get3DIntersectionPoint(CanvasPoint(x, y, focalLength), cameraPosition, focalLength, cameraOrientation) - cameraPosition);
			RayTriangleIntersection triangleInter = getClosestIntersection(cameraPosition, rayDirection, v, noIntTriangles, barycentric);
			ModelTriangle interTriangle = triangleInter.intersectedTriangle;
			if (triangleInter.triangleIndex == -1) continue; // nothing to draw


			std::pair<glm::vec3, glm::vec2> NORMALES = getNormal(barycentric, vertexNormals, trToNormal, triangleInter.intersectedTriangle, v);
			glm::vec3 NORMAL = NORMALES.first;
			glm::vec3 NORMAL2 = NORMAL;
			glm::vec2 texture = NORMALES.second;
			
			//reflection
			if (triangleInter.intersectedTriangle.colour.name == "Blue" && REFLEC) {
				glm::vec3 Ri = rayDirection;
				glm::vec3 reflection = Ri - 2.0f * NORMAL * (glm::dot(Ri, NORMAL));
				int depth = 0;
				std::pair<Colour, std::pair<std::array<TexturePoint, 3>, glm::vec2> >miRRor = mirror(triangleInter.intersectionPoint, reflection, v, depth, vertexNormals, trToNormal, NORMAL, triangleInter.intersectedTriangle);
				triangleInter.intersectedTriangle = interTriangle;
				triangleInter.intersectedTriangle.colour = miRRor.first;
				texture = miRRor.second.second;
				triangleInter.intersectedTriangle.texturePoints = miRRor.second.first;
				NORMAL = NORMAL2;
			}
			
			if (triangleInter.intersectedTriangle.colour.name == "Red" && REFLEC2) {
				
				glm::vec3 refraction = refract(rayDirection, NORMAL, 1.5f);
				int depth = 0;
				std::pair<Colour, std::pair<std::array<TexturePoint, 3>, glm::vec2> >miRRor = refra(triangleInter.intersectionPoint, refraction, v, depth, vertexNormals, trToNormal, NORMAL, triangleInter.intersectedTriangle, true);
				triangleInter.intersectedTriangle = interTriangle;
				triangleInter.intersectedTriangle.colour = miRRor.first;
				texture = miRRor.second.second;
				triangleInter.intersectedTriangle.texturePoints = miRRor.second.first;
				NORMAL = NORMAL2;
			}
			
			glm::vec3 lightDirection = glm::normalize(lightLocation -  triangleInter.intersectionPoint);
			noIntTriangles = 0;
			RayTriangleIntersection lightInter = getClosestIntersection(triangleInter.intersectionPoint + triangleInter.intersectedTriangle.normal * 0.0001f, lightDirection, v, noIntTriangles, empty);
			Colour triangleColour;
			
			float all_speculars = 0;
			float all_diffuse = 0;
			float r;
			for (auto &light: lights) {
				lightDirection = glm::normalize(light -  triangleInter.intersectionPoint);
				noIntTriangles = 0;
				r =  glm::length(light - triangleInter.intersectionPoint);
				lightInter = getClosestIntersection(triangleInter.intersectionPoint + triangleInter.intersectedTriangle.normal * 0.0001f, lightDirection, v, noIntTriangles, empty);
				
				if (lightInter.distanceFromCamera > r) {
					all_speculars += getSpecular(lightDirection, NORMAL, cameraPosition, triangleInter.intersectionPoint, 8);
					all_diffuse += getDiffuse(r, NORMAL, lightDirection);
				}

			}

			//if textureON
			if (triangleInter.intersectedTriangle.texturePoints[0].x != 0 || triangleInter.intersectedTriangle.texturePoints[0].y != 0 || triangleInter.intersectedTriangle.texturePoints[1].x != 0 || triangleInter.intersectedTriangle.texturePoints[1].y != 0 || triangleInter.intersectedTriangle.texturePoints[2].x != 0 || triangleInter.intersectedTriangle.texturePoints[2].y != 0) {
				uint32_t c = MAP.pixels[int(texture.x) + int(texture.y) * MAP.width];
				int r = (c >> 16) & 0xff; // red
 				int g = (c >> 8) & 0xff; // green
 				int b = c  & 0xff; // blue
				triangleColour = Colour(r, g, b);
			} 
			else {
				triangleColour = triangleInter.intersectedTriangle.colour;
			}

			//sphere.obj
			if (SPHERE)triangleColour = Colour(255, 0, 0);

			float diffuse = all_diffuse / lights.size();
			float specular = all_speculars / lights.size();
			float ambientStrength = 0.5;
	
			
			window.setPixelColour(x, y, getColour(triangleColour.red, triangleColour.green, triangleColour.blue, diffuse, specular, ambientStrength));
		}
	} 

}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &cameraPosition, glm::mat3 &cameraOrientation, std::vector <ModelTriangle> &v, float focalLength, glm::vec3 &lightLocation, std::unordered_map <std::string, TextureMap> textures, std::vector <glm::vec3> &vertexNormals, std::vector < std::vector <int> > &trToNormal, std::vector <glm::vec3> &lights, bool &animation) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "LEFT" << std::endl;
			cameraPosition[0] -= 0.1;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			std::cout << "RIGHT" << std::endl;
			cameraPosition[0] += 0.1;
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			std::cout << "UP" << std::endl;
			cameraPosition[1] += 0.1;
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			std::cout << "DOWN" << std::endl;
			cameraPosition[1] -= 0.1;
		}
		else if (event.key.keysym.sym == SDLK_z) {
			std::cout << "BACKWARD" << std::endl;
			cameraPosition[2] += 0.1;
		}
		else if (event.key.keysym.sym == SDLK_x) {
			std::cout << "FORWARD" << std::endl;
			cameraPosition[2] -= 0.1;
		}
		else if (event.key.keysym.sym == SDLK_s) { 
			std::cout << "Change camera position" << std::endl;
			float rad = 0.10;
			glm::mat3 xRotation = glm::mat3(
				1, 0, 0,
				0, cos(rad), sin(rad),
				0, (-1) * sin(rad), cos(rad)
			);

			cameraPosition = cameraPosition * xRotation;

		}
		else if (event.key.keysym.sym == SDLK_w) { 
			std::cout << "Change camera position" << std::endl;
			float rad = -0.10;
			glm::mat3 xRotation = glm::mat3(
				1, 0, 0,
				0, cos(rad), sin(rad),
				0, (-1) * sin(rad), cos(rad)
			);

			cameraPosition = cameraPosition * xRotation;
		}
		else if (event.key.keysym.sym == SDLK_d) { 
			std::cout << "Change camera position" << std::endl;
			float rad = 0.10;
			glm::mat3 yRotation = glm::mat3(
				cos(rad), 0, -sin(rad),
				0, 1, 0,
				sin(rad), 0, cos(rad)
			);

			cameraPosition = cameraPosition * yRotation;
		}
		else if (event.key.keysym.sym == SDLK_a) { 
			std::cout << "Change camera position" << std::endl;
			float rad = -0.10;
			glm::mat3 yRotation = glm::mat3(
				cos(rad), 0, -sin(rad),
				0, 1, 0,
				sin(rad), 0, cos(rad)
			);

			cameraPosition = cameraPosition * yRotation;
		}
		else if (event.key.keysym.sym == SDLK_k) { 
			std::cout << "Change camera orientation" << std::endl;
			float rad = 0.10;
			glm::mat3 xRotation = glm::mat3(
				1, 0, 0,
				0, cos(rad), sin(rad),
				0, (-1) * sin(rad), cos(rad)
			);

			cameraOrientation = cameraOrientation * xRotation;
		}
		else if (event.key.keysym.sym == SDLK_i) { 
			std::cout << "Change camera orientation" << std::endl;
			float rad = -0.10;
			glm::mat3 xRotation = glm::mat3(
				1, 0, 0,
				0, cos(rad), sin(rad),
				0, (-1) * sin(rad), cos(rad)
			);

			cameraOrientation = cameraOrientation * xRotation;
		}
		else if (event.key.keysym.sym == SDLK_l) { 
			std::cout << "Change camera orientation" << std::endl;
			float rad = 0.10;
			glm::mat3 yRotation = glm::mat3(
				cos(rad), 0, -sin(rad),
				0, 1, 0,
				sin(rad), 0, cos(rad)
			);

			cameraOrientation = cameraOrientation * yRotation;
		}
		else if (event.key.keysym.sym == SDLK_j) { 
			std::cout << "Change camera orientation" << std::endl;
			float rad = -0.10;
			glm::mat3 yRotation = glm::mat3(
				cos(rad), 0, -sin(rad),
				0, 1, 0,
				sin(rad), 0, cos(rad)
			);

			cameraOrientation = cameraOrientation * yRotation;
		}

		else if (event.key.keysym.sym == SDLK_1) {
			std::cout << "Render wireframed scene" << std::endl;
			std::vector <std::vector <float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0.0)); //not needed for wireframed - remains default
			float speed = 0.01;
			drawWireframeScene(window, v, cameraPosition, focalLength, cameraOrientation, depthBuffer, speed, false);

		}
		else if (event.key.keysym.sym == SDLK_2) {
			std::cout << "Render rasterised scene" << std::endl;
			std::vector <std::vector <float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0.0));	
			float speed = 0.01;
			drawRasterisedScene(window, v, cameraPosition, focalLength, depthBuffer, cameraOrientation, textures, speed, false);
		}
		else if (event.key.keysym.sym == SDLK_3) {
			std::cout << "Render raytraced scene" << std::endl;
			drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, v, lightLocation, vertexNormals, trToNormal, textures["Cobbles"], lights,  false, false, false, false,false);
		}
		else if (event.key.keysym.sym == SDLK_p) {
			std::cout<< "Move source light1" << std::endl;
			lightLocation[1] += 0.1;
		}
		else if (event.key.keysym.sym == SDLK_o) {
			std::cout<< "Move source light2" << std::endl;
			lightLocation[1] -= 0.1;
			std::cout << lightLocation[1] << " ";
		}
		else if (event.key.keysym.sym == SDLK_0) {
			std::cout<< "Start Animation" << std::endl;
			animation = true;
			cameraPosition = glm::vec3(0.0, 0.0, 3.0);
			cameraOrientation =glm::mat3(
			1, 0, 0, //RIGHT
			0, 1, 0, //UP
			0, 0, 1 //FORWARD
			);
		}
	} else if (event.type == SDLK_4) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	glm::vec3 cameraPositionSphere(0.0, 30.0, 15); //spehere.obj
	glm::vec3 cameraPosition(0.0, 0.0, 3.0);  //++cornell
	glm::vec3 initialCamPos = cameraPosition;

	glm::mat3 cameraOrientation = glm::mat3(
			1, 0, 0, //RIGHT
			0, 1, 0, //UP
			0, 0, 1 //FORWARD
		);
	glm::mat3 initialCamOr = cameraOrientation;
	std::unordered_map <std::string, TextureMap> textures;
	std::vector < std::vector<int> > trToNormal;
	std::vector < glm::vec3 > vertexNormals;
	TextureMap MAP = TextureMap("texture.ppm");
	TextureMap MAP2 = TextureMap("check3.ppm");
	float focalLength = 2.0;

	std::vector <ModelTriangle> arr = parseObj("textured-cornell-box.obj", 0.17, textures, vertexNormals, trToNormal, MAP);
	//std::vector <ModelTriangle> arr = parseObj("cornell-box.obj", 0.17, textures, vertexNormals, trToNormal, MAP);
	//std::vector <ModelTriangle> arr = parseObj("sphere2.obj", 0.17, textures, vertexNormals, trToNormal, MAP);

	glm::vec3 lightLocation2(0.0, 1.3, 1.2); //sphere.obj
	glm::vec3 lightLocation(0.0, 0.20, 0.20); //textured cornell
	std::vector <glm::vec3> lights, lightsSph;
	lightsSph.push_back(lightLocation2);
	lights.push_back(lightLocation);
	float xx = 0.0, yy = 0.20, zz = 0.2, t = 1;
	for (int i = 0; i < 2; ++i) {
		lights.push_back(glm::vec3(xx - 0.20 * t, yy, zz));
		lights.push_back(glm::vec3(xx + 0.20 * t, yy, zz));
		lights.push_back(glm::vec3(xx, yy - 0.2 * t, zz));
		lights.push_back(glm::vec3(xx, yy + 0.2 * t, zz));
		t++;
	}

	bool animation = false;
	bool wf = true, rst = false, rt = false, rt2 = false, rt3 = false, sph = false, zoom = false;
	float speed1 = 0.01, speed2 = 0.1;
	int nr = 2 * 3.14 / speed1, halfNr = nr / 2;
	//press 0 to start animation (one time)
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window, cameraPosition, cameraOrientation, arr, focalLength, lightLocation, textures, vertexNormals, trToNormal ,lights, animation);
		std::vector <std::vector <float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0.0));	
		if (!animation) {
			drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP2, lights, true, false, zoom, false, false);
		} 
		if (wf && nr && animation) {
			drawWireframeScene(window, arr, cameraPosition, focalLength, cameraOrientation, depthBuffer, speed1, zoom);
			//drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP, lights);
			nr--;
			if (nr <= halfNr) zoom = true;

			if (nr == 0) {
				rst = true;
				wf = false;
				zoom = false;
				nr = 2 * 3.14 / speed1;
				halfNr = nr / 2;
				cameraPosition = initialCamPos;
				cameraOrientation = initialCamOr;
			}
		}
		else
		if (rst && nr) {
			drawRasterisedScene(window, arr, cameraPosition, focalLength, depthBuffer, cameraOrientation, textures, speed1, zoom);
			nr--;
			if (nr <= halfNr) zoom = true;

			if (nr == 0) {
				rst = false;
				rt = true;
				zoom = true;
				nr =(2 * 3.14 / speed2) + 1;
				halfNr = nr / 2;
				cameraPosition = initialCamPos;
				cameraOrientation = initialCamOr;
			}
		}
		else
		if (rt && nr) {
			drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP2, lights, false, false, zoom, sph, true);
			nr--;
			if (nr <= halfNr) zoom = false;

			if (nr == 0) {
				rt = false;
				rt2 = true;
				zoom = true;
				nr =(2 * 3.14 / speed2) + 1;
				halfNr = nr / 2;
				cameraPosition = initialCamPos;
				cameraOrientation = initialCamOr;
			}
		}
		else if (rt2 && nr) {
			drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP2, lights, true, false, zoom, sph, true);
			nr--;
			if (nr <= halfNr) zoom = false;
			if (nr == 0) {
				rt2 = false;
				rt3 = true;
				zoom = true;
				nr =(2 * 3.14 / speed2) + 1;
				halfNr = nr / 2;
				cameraPosition = initialCamPos;
				cameraOrientation = initialCamOr;
			}
		}
		else if (rt3 && nr) {
			drawRayTracedScene(window, cameraPosition, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP2, lights, false, true, zoom, sph, true);
			nr--;
			if (nr <= halfNr) zoom = false;
			if (nr == 0) {
				rt3 = false;
				sph = true;
				nr =(2 * 3.14 / speed2) + 1;
				cameraOrientation = initialCamOr;
				vertexNormals.clear();
				trToNormal.clear();
				arr = parseObj("sphere2.obj", 0.17, textures, vertexNormals, trToNormal, MAP2);
			}
		}
		else if(sph && nr) {
			drawRayTracedScene(window, cameraPositionSphere, focalLength, cameraOrientation, arr, lightLocation, vertexNormals, trToNormal, MAP, lightsSph, false, false, zoom, sph, true);
			nr--;
		}


		window.renderFrame();
	}


}
