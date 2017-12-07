#include <pathgen/PathGen.h>
#include <pangolin/pangolin.h>
#include <SceneGraph/SceneGraph.h>
#include <SceneGraph/GLDynamicGrid.h>
#include <glog/logging.h>
#include <gflags/gflags.h>

/*-----------COMMAND LINE FLAGS-----------------------------------------------*/
DEFINE_bool(save_imu_files, false, "save deubg imu measurement, poses etc files");
DEFINE_bool(save_cam_files, false, "save deubg cam poses");
DEFINE_bool(display, true, "display poses with pangolin");
DEFINE_uint32(imu_rate, 100, "imu frequency");
DEFINE_uint32(camera_rate, 5, "camera frequency");
DEFINE_double(duration, 10.0, "duration of motion");

/*----------------------------------------------------------------------------*/

struct GuiVars {
    pangolin::OpenGlRenderState render_state;
    int image_width;
    int image_height;
    SceneGraph::GLSceneGraph scene_graph;
    SceneGraph::GLDynamicGrid grid;
    SceneGraph::GLText text;
    pangolin::View* grid_view;
    pangolin::OpenGlRenderState gl_render3d;
    std::unique_ptr<SceneGraph::HandlerSceneGraph> sg_handler_;
    pangolin::View* panel_view;
};

GuiVars gui_vars_;
std::vector<pangolin::GlTexture> gl_tex_;
SceneGraph::AxisAlignedBoundingBox aabb_;
std::vector<std::unique_ptr<SceneGraph::GLAxis> > axes_;
std::shared_ptr<SceneGraph::GLPrimitives<>> line_strip_;

std::unique_ptr<pangolin::Var<bool> > random_seed_;
std::unique_ptr<pangolin::Var<bool> > plot_imu_data_;
std::unique_ptr<pangolin::Var<bool> > run_;
std::unique_ptr<pangolin::Var<int> > distance_;
std::unique_ptr<pangolin::Var<int> > radius_;


void InitGui()
{

    const int window_width = 1024;
    const int window_height = 800;

    pangolin::CreateGlutWindowAndBind("PathGen", window_width,
                                      window_height);

    gui_vars_.render_state.SetModelViewMatrix(pangolin::IdentityMatrix());
    gui_vars_.render_state.SetProjectionMatrix(
                pangolin::ProjectionMatrixOrthographic(
                    0, window_width, 0, window_height, 0, 1000));

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    gui_vars_.scene_graph.AddChild(&gui_vars_.grid);

    gui_vars_.scene_graph.AddChild(&gui_vars_.text);


    gui_vars_.grid_view = &pangolin::Display("grid")
            .SetAspect(-(float)window_width / (float)window_height);

    gui_vars_.panel_view = &pangolin::CreatePanel("ui").
            SetBounds(0, 1.0, 0, pangolin::Attach::Pix(100));
    random_seed_.reset(new pangolin::Var<bool>("ui.RandomSeed", false,true));
    plot_imu_data_.reset(new pangolin::Var<bool>("ui.RandomSeed", true,true));
    run_.reset(new pangolin::Var<bool>("ui.Run", false,false));
    distance_.reset(new pangolin::Var<int>("ui.Distance", 10, 2, 20));
    radius_.reset(new pangolin::Var<int>("ui.Radius", 10, 2, 20));

    gui_vars_.gl_render3d.SetProjectionMatrix(
                pangolin::ProjectionMatrix(window_width, window_height,
                                           420, 420, 320, 240, 0.01, 5000));

    gui_vars_.gl_render3d.SetModelViewMatrix(
                pangolin::ModelViewLookAt(3, 3, 4, 0, 0, 0,
                                          pangolin::AxisZ));

    gui_vars_.sg_handler_.reset(new SceneGraph::HandlerSceneGraph(
                                    gui_vars_.scene_graph,
                                    gui_vars_.gl_render3d,
                                    pangolin::AxisZ,
                                    50.0f));

    gui_vars_.grid_view->SetHandler(gui_vars_.sg_handler_.get());
    gui_vars_.grid_view->SetDrawFunction(
                SceneGraph::ActivateDrawFunctor(
                    gui_vars_.scene_graph,
                    gui_vars_.gl_render3d));

    pangolin::Display("multi").SetBounds(0, 1.0, pangolin::Attach::Pix(100), 1.0).SetLayout(
                pangolin::LayoutEqual);

    pangolin::Display("multi").AddDisplay(*gui_vars_.grid_view);

    SceneGraph::GLSceneGraph::ApplyPreferredGlSettings();
    glClearColor(0.0, 0.0, 0.0, 1.0);

    line_strip_.reset(new SceneGraph::GLPrimitives<>);

    gui_vars_.scene_graph.AddChild(line_strip_.get());
}

using namespace pathgen;

std::shared_ptr<PoseSpline> Run(std::vector<PosePtr>* poses, std::vector<PosePtr>* imu_poses,
                                ImuMeasurementDeque* imu_measurements){

//    uint64_t seed = 0.0;
//    if(*random_seed_){
//        seed = (uint64_t)time(0);
//    }
//    srand(seed);

//    // generate a discrete set of poses
//    *poses = PathGenerator::RandomWalk(10, *radius_, *distance_);

//    // now run a continous spline through the poses
//    std::shared_ptr<PoseSpline> spline = PathGenerator::GenerateSplineFromPoses(*poses);

//    // generate imu measurements corresponding to the spline
//    ImuParameters imu_parameters;
//    ImuPathOptions path_options;
//    path_options.add_noise_to_imu_measurements = true;
//    path_options.duration = 10;
//    IMUGenerator::GenerateInertialMeasurements(*spline, imu_parameters,
//                                               path_options,
//                                               imu_measurements,
//                                               imu_poses);

//    return spline;
}

//bool ComputeJacobian(const double* x,
//                                                 double* jacobian) const {
//  jacobian[0] = -x[1]; jacobian[1]  = -x[2]; jacobian[2]  = -x[3];  // NOLINT
//  jacobian[3] =  x[0]; jacobian[4]  =  x[3]; jacobian[5]  = -x[2];  // NOLINT
//  jacobian[6] = -x[3]; jacobian[7]  =  x[0]; jacobian[8]  =  x[1];  // NOLINT
//  jacobian[9] =  x[2]; jacobian[10] = -x[1]; jacobian[11] =  x[0];  // NOLINT
//  return true;
//}

bool ComputeJacobian(const double* x,double* jacobian) {
  jacobian[0] =  x[3]; jacobian[1]  =  x[2]; jacobian[2]  = -x[1];  // NOLINT
  jacobian[3] = -x[2]; jacobian[4]  =  x[3]; jacobian[5]  =  x[0];  // NOLINT
  jacobian[6] =  x[1]; jacobian[7]  = -x[0]; jacobian[8]  =  x[3];  // NOLINT
  jacobian[9] = -x[0]; jacobian[10] = -x[1]; jacobian[11] = -x[2];  // NOLINT
  return true;
}

int main(int argc, char** argv)
{
    google::InitGoogleLogging(argv[0]);
    google::ParseCommandLineFlags(&argc, &argv, true);


    // set up logging to console
    FLAGS_stderrthreshold = 0;  // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;
    FLAGS_logtostderr = 1;

    InitGui();

    uint64_t seed = (uint64_t)time(0);
    srand(seed);


    // generate a discrete set of poses
    std::vector<PosePtr> poses;
    ImuMeasurementDeque imu_measurements;
    std::vector<PosePtr> imu_poses;

    std::shared_ptr<PoseSpline> spline =
            Run(&poses, &imu_poses, &imu_measurements);
    // now run a continous spline through the poses

    // generate imu measurements corresponding to the spline

    if(FLAGS_display){

        while(!pangolin::ShouldQuit()){
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glColor4f(1.0f,1.0f,1.0f,1.0f);

            if (pangolin::Pushed(*run_)) {
                poses.clear();
                imu_poses.clear();
                imu_measurements.clear();
                for(const auto& a : axes_){
                    gui_vars_.scene_graph.RemoveChild(a.get());
                }
                axes_.clear();
                gui_vars_.scene_graph.RemoveChild(line_strip_.get());
                line_strip_->Clear();
                Run(&poses, &imu_poses, &imu_measurements);
            }


            // activate the 3d view
            gui_vars_.grid_view->ActivateAndScissor(gui_vars_.gl_render3d);

            if(axes_.size() == 0){
                for(const PosePtr& pose : imu_poses){
                    axes_.push_back(std::unique_ptr<SceneGraph::GLAxis>(
                                        new SceneGraph::GLAxis(0.05)));
                    axes_.back()->SetPose(pose->matrix());
                    aabb_.Insert(pose->translation());
                    gui_vars_.grid.set_bounds(aabb_);
                    gui_vars_.scene_graph.AddChild(axes_.back().get());

                    // add a line connecting the frames
                    Eigen::Vector3f vertex = pose->translation().cast<float>();
                    line_strip_->AddVertex(vertex);
                }
            }

            pangolin::FinishFrame();
        }

    }

    if(FLAGS_save_imu_files){
        // save the imu measurements to a file
        IMUGenerator::SaveMeasurementFiles(imu_measurements);

        // also save the imu poses to a file
        SavePoses(imu_poses, "imu_poses.txt");
    }

    if(FLAGS_save_cam_files){
        SaveCameraPoses(*spline, FLAGS_camera_rate,
                        imu_measurements.front().timestamp);
    }



}
