//#include "phasor.h"
//
//
//std::vector<double> phasor::compute_phasor(
//        unsigned short *micro_times,
//        std::vector<int> &idxs,
//        double frequency,
//        int minimum_number_of_photons,
//        double g_irf,
//        double s_irf
//){
//    double factor = (2. * frequency * M_PI);
//    std::vector<double> re{-1, -1};
//    double g_sum = 0.0;
//    double s_sum = 0.0;
//    double sum = (double) idxs.size();
//    for(auto &idx: idxs){
//        auto mt = micro_times[idx];
//        g_sum += std::cos(mt * factor);
//        s_sum += std::sin(mt * factor);
//    }
//    if(sum > minimum_number_of_photons){
//        double g_exp = g_sum / std::max(1., sum);
//        double s_exp = s_sum / std::max(1., sum);
//        re[0] = phasor::g(g_irf, s_irf, g_exp, s_exp);
//        re[1] = phasor::s(g_irf, s_irf, g_exp, s_exp);
//    }
//    return re;
//}
//
//
//std::vector<double> phasor::compute_phasor_all(
//        unsigned short* microtimes, int n_microtimes,
//        double frequency
//){
//    double factor = (2. * frequency * M_PI);
//    double sum = n_microtimes;
//    double g_sum = 0.0; double s_sum = 0.0;
//    for(int i=0;i<n_microtimes;i++){
//        g_sum += std::cos(microtimes[i] * factor);
//        s_sum += std::sin(microtimes[i] * factor);
//    }
//    std::vector<double> re{g_sum / sum, s_sum / sum};
//    return re;
//}
//
//
//double phasor::g(
//            double g_irf, double s_irf,
//            double g_exp, double s_exp
//    ) {
//        return 1. / (g_irf * g_irf + s_irf * s_irf) * (g_irf * g_exp + s_irf * s_exp);
//    }
//
//    double phasor::s(
//            double g_irf, double s_irf,
//            double g_exp, double s_exp
//    ) {
//        return 1. / (g_irf * g_irf + s_irf * s_irf) * (g_irf * s_exp - s_irf * g_exp);
//    }
//
//
//void phasor::get_phasor_image(
//        float **output, int *dim1, int *dim2, int *dim3, int *dim4,
//        std::shared_ptr<CLSMImage> image,
//        TTTR *tttr_irf,
//        double frequency,
//        int minimum_number_of_photons,
//        bool stack_frames
//) {
//    double g_irf=1.0, s_irf=0.0;
//    // compute g, s phasor of IRF
//    if(tttr_irf!= nullptr){
//        unsigned short *micro_times; int n_micro_times;
//        tttr_irf->get_micro_time(&micro_times, &n_micro_times);
//        std::vector<double> gs = phasor::compute_phasor_all(
//                micro_times, n_micro_times,
//                frequency);
//        g_irf = gs[0];
//        s_irf = gs[1];
//    }
//
//    auto tttr_data = image->tttr;
//    // if stack frames only one output frame
//    int o_frames = stack_frames? 1: image->get_n_frames();
//    auto header = tttr_data->get_header();
//    if(frequency<0){
//        frequency = 1. / (1000 / header->get_micro_time_resolution());
//    }
//    double factor = (2. * frequency * M_PI);
//
//#if VERBOSE_FIT2X
//    std::clog << "GET_PHASOR_IMAGE..." << std::endl;
//    std::clog << "-- frequency [GHz]: " << frequency << std::endl;
//    std::clog << "-- stack_frames: " << stack_frames << std::endl;
//    std::clog << "-- minimum_number_of_photons: " << minimum_number_of_photons << std::endl;
//#endif
//
//    int n_frames = image->get_n_frames();
//    int n_pixel = image->get_n_pixel();
//    int n_lines = image->get_n_lines();
//    //
//    unsigned short *micro_times; int n_micro_times;
//    tttr_data->get_micro_time(&micro_times, &n_micro_times);
//
//    // allocate the output array
//    auto* t = (float *) calloc(o_frames * n_lines * n_pixel * 2, sizeof(float));
//    size_t pixel_nbr = 0;
//    for(int i_line=0; i_line < n_lines; i_line++){
//        for(int i_pixel=0; i_pixel < n_pixel; i_pixel++){
//            if(stack_frames){
//                std::vector<int> idxs = {};
//                pixel_nbr = i_line  * (image->get_n_pixel() * 2) + i_pixel * 2;
//                for(int i_frame=0; i_frame < n_frames; i_frame++){
//                    auto frame = image->get_frames()[i_frame];
//                    auto line = frame->get_lines()[i_line];
//                    auto idx = line->get_pixels()[i_pixel]->_tttr_indices;
//                    idxs.insert(idxs.end(), idx.begin(), idx.end());
//                }
//                auto r = phasor::compute_phasor(micro_times, idxs, frequency, minimum_number_of_photons, g_irf, s_irf);
//                t[pixel_nbr + 0] = r[0];
//                t[pixel_nbr + 1] = r[1];
//            } else{
//                for(int i_frame=0; i_frame < n_frames; i_frame++){
//                    pixel_nbr = i_frame * (n_lines * n_pixel * 2) + i_line  * (n_pixel * 2) + i_pixel * 2;
//                    auto frame = image->get_frames()[i_frame];
//                    auto line = frame->get_lines()[i_line];
//                    auto idx = line->get_pixels()[i_pixel]->_tttr_indices;
//                    auto r = phasor::compute_phasor(micro_times, idx, frequency, minimum_number_of_photons, g_irf, s_irf);
//                    t[pixel_nbr + 0] = r[0];
//                    t[pixel_nbr + 1] = r[1];
//                }
//            }
//        }
//    }
//    *dim1 = (int) o_frames;
//    *dim2 = (int) n_lines;
//    *dim3 = (int) n_pixel;
//    *dim4 = (int) 2;
//    *output = t;
//}