#include "config.hpp"
#include "field-math.hpp"
#include "optimizer.hpp"
#include "parametrizer.hpp"
#include <stdlib.h>

#ifdef WITH_CUDA
#include <cuda_runtime.h>
#endif

using namespace qflow;

Parametrizer field;

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

#ifdef WITH_CUDA
    cudaFree(0);
#endif
    int t1, t2;
    std::string input_obj, output_obj;
    int faces = -1;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-f") == 0) {
            sscanf(argv[i + 1], "%d", &faces);
        } else if (strcmp(argv[i], "-i") == 0) {
            input_obj = argv[i + 1];
        } else if (strcmp(argv[i], "-o") == 0) {
            output_obj = argv[i + 1];
        } else if (strcmp(argv[i], "-sharp") == 0) {
            field.flag_preserve_sharp = 1;
        } else if (strcmp(argv[i], "-adaptive") == 0) {
            field.flag_adaptive_scale = 1;
        } else if (strcmp(argv[i], "-mcf") == 0) {
            field.flag_minimum_cost_flow = 1;
        } else if (strcmp(argv[i], "-sat") == 0) {
            field.flag_aggresive_sat = 1;
        } else if (strcmp(argv[i], "-seed") == 0) {
            field.hierarchy.rng_seed = atoi(argv[i + 1]);
        }
    }
    printf("%d %s %s\n", faces, input_obj.c_str(), output_obj.c_str());
    if (input_obj.size() >= 1) {
        field.Load(input_obj.c_str());
    } else {
        assert(0);
        // field.Load((std::string(DATA_PATH) + "/fertility.obj").c_str());
    }

    printf("Initialize...\n");
    t1 = GetCurrentTime64();
    field.Initialize(faces);
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    printf("Solve Orientation Field...\n");
    t1 = GetCurrentTime64();

    Optimizer::optimize_orientations(field.hierarchy);
    field.ComputeOrientationSingularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    if (field.flag_adaptive_scale == 1) {
        printf("Estimate Slop...\n");
        t1 = GetCurrentTime64();
        field.EstimateSlope();
        t2 = GetCurrentTime64();
        printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    }
    printf("Solve for scale...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_scale(field.hierarchy, field.rho, field.flag_adaptive_scale);
    field.flag_adaptive_scale = 1;
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    printf("Solve for position field...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_positions(field.hierarchy, field.flag_adaptive_scale);

    field.ComputePositionSingularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    t1 = GetCurrentTime64();
    printf("Solve index map...\n");
    field.ComputeIndexMap();
    t2 = GetCurrentTime64();
    printf("Indexmap Use %lf seconds\n", (t2 - t1) * 1e-3);
    printf("Writing the file...\n");

    if (output_obj.size() < 1) {
        assert(0);
        // field.OutputMesh((std::string(DATA_PATH) + "/result.obj").c_str());
    } else {
        field.OutputMesh(output_obj.c_str());
    }
    printf("finish...\n");
    //	field.LoopFace(2);
    return 0;
}
