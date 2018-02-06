#if defined __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "config.hpp"
#include "parametrizer.hpp"
#include "optimizer.hpp"
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include "gldraw.hpp"
#include "field-math.hpp"


#ifdef WITH_CUDA
#include <cuda_runtime.h>
#endif

Parametrizer field;

int main(int argc, char** argv)
{
#ifdef WITH_SCALE
	int with_scale = 1;
#else
	int with_scale = 0;
#endif

#ifdef WITH_CUDA
	cudaFree(0);
#endif
	int t1, t2;
	std::string input_obj, output_obj;
	int faces = -1;
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], "-f")==0) {
			sscanf(argv[i + 1], "%d", &faces);
		}
		if (strcmp(argv[i], "-i") == 0) {
			input_obj = argv[i + 1];
		}
		if (strcmp(argv[i], "-o") == 0) {
			output_obj = argv[i + 1];
		}
	}
	printf("%d %s %s\n", faces, input_obj.c_str(), output_obj.c_str());
	fflush(stdout);
	if (input_obj.size() >= 1)
        field.Load(input_obj.c_str());
    else
        field.Load((std::string(DATA_PATH) + "/fertility.obj").c_str());
	
	printf("Initialize...\n");
	fflush(stdout);
	field.Initialize(faces, with_scale);

	printf("Solve Orientation Field...\n");
	fflush(stdout);
	t1 = GetCurrentTime64();
	
	Optimizer::optimize_orientations(field.hierarchy);
	field.ComputeOrientationSingularities();
	t2 = GetCurrentTime64();
	printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
	fflush(stdout);

	if (with_scale == 1) {
		printf("estimate for scale...\n");
		fflush(stdout);
		t1 = GetCurrentTime64();
		field.EstimateScale();
		t2 = GetCurrentTime64();
		printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

		printf("Solve for scale...\n");
		t1 = GetCurrentTime64();
		Optimizer::optimize_scale(field.hierarchy);
		t2 = GetCurrentTime64();
		printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
	}

	printf("Solve for position field...\n");
	fflush(stdout);
	t1 = GetCurrentTime64();
	Optimizer::optimize_positions(field.hierarchy, with_scale);
	field.ComputePositionSingularities(with_scale);
	t2 = GetCurrentTime64();
	printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
	fflush(stdout);
	t1 = GetCurrentTime64();
	printf("Solve index map...\n");
	fflush(stdout);
	field.ComputeIndexMap(with_scale);
	t2 = GetCurrentTime64();
	printf("Indexmap Use %lf seconds\n", (t2 - t1) * 1e-3);
	fflush(stdout);

	printf("Writing the file...\n");
	fflush(stdout);
	if (output_obj.size() < 1)
		field.ExtractMesh((std::string(DATA_PATH) + "/result.obj").c_str());
	else
		field.ExtractMesh(output_obj.c_str());
	printf("finish...\n");
	fflush(stdout);
	//	field.LoopFace(2);
#ifdef WITH_OPENGL
	gldraw();
#endif
	return 0;
}
