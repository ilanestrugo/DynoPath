//
//  main.cpp
//  Sorting code (Simulation and algorithm)
//
//  Created by Tal Raviv on 19/02/2024.
//


#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <chrono>
#include <limits>
#include <list>
#include <tuple>
#include <tuple>
#include <random>
#include "argparse.hpp"

//init 5 integers for 4 possible directions + 1 for stay in place
const int dir_stay = 0;
const int dir_up = 1;
const int dir_right = 2;
const int dir_down = 3;
const int dir_left = 4;

//define TEG vector
std::vector<char> TEG; 
// We could use bool to save memory but it will take 20% more time to update
// We keep this as a global variable to save the need § pass it many times to find_path, get_TEG, and set_TEG

//define global Parameters
int Lx;
int Ly;
int T;   // number of takts in the simulation not including warmup and cool down periods
int T1;  // number of takts in the simulation including warmup and cool down periods
int cool_down;
bool verbose = false;
bool do_animation = false;

int count_entered = 0;
std::vector<std::list<std::tuple<int, std::tuple<int,int>, std::tuple<int,int>, int>>> AnimationData;

argparse::ArgumentParser program("sorting");

//function to get data from the 4-dim array that represents TEG
bool get_TEG(int x, int y,int t, int d) {
    
#ifdef QA
    return TEG.at((((x * Ly + y) * T1) + t)*5 + d);
#endif

#ifndef QA
    return TEG[(((x * Ly + y) * T1) + t)*5 + d];
#endif
    
}

//function to update the 4-dim array that represents TEG
void set_TEG(int x, int y,int t, int d, bool v) {
#ifdef QA
    if (x>= Lx or x< 0 or y < 0 or y>= Ly or t<0 or t>=T1 or d<0 or d>4) {
        std::cout << "Panic: out of range\n";
        exit(1);
    }
        
    TEG.at((((x * Ly + y) * T1) + t)*5 + d) = v;
#endif

#ifndef QA
    TEG[(((x * Ly + y) * T1) + t)*5 + d] = v;
#endif

}

// function to faciltate acess to the 2d and 3d arrays
int c3d(int x,int y,int t) {
    return ((x*Ly)+y)*cool_down+t;
}

int c2d(int x,int y) {
    return (x*Ly)+y;
}

//Function to split a string by a substring
std::vector<std::string> splitString(const std::string& input, const std::string& delimiter) {
    std::vector<std::string> tokens;

    size_t start = 0, end = 0;
    while ((end = input.find(delimiter, start)) != std::string::npos) {
        tokens.push_back(input.substr(start, end - start));
        start = end + delimiter.length();
    }

    // Add the last token
    tokens.push_back(input.substr(start));

    return tokens;
}

//functions used for statistics
double calculateMean(const std::vector<double_t>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}

double calculateStandardDeviation(const std::vector<double_t>& data, double &mean) {
    // Calculate the mean
    mean = calculateMean(data);

    // Calculate the squared differences
    double squaredDifferencesSum = 0.0;
    for (double value : data) {
        double difference = value - mean;
        squaredDifferencesSum += difference * difference;
    }

    // Calculate the mean of squared differences
    double meanOfSquaredDifferences = squaredDifferencesSum / data.size();

    // Take the square root to get the standard deviation
    double standardDeviation = std::sqrt(meanOfSquaredDifferences);

    return standardDeviation;
}

// this function find the shortest path on the TEG and update the TEG acordignly
// the function returns the number of moves require to acomplishe the retrival
// the flow time value is returned in a by referance argument
int find_path(int t_base, int x0, int y0, int dest, std::list<std::tuple<int, int, int>> OutputsWithDest, int &flow_time, long int &_total_nodes) {
    std::vector<char> from_node(Lx*Ly*cool_down); // value 0-4 - direction from which we arrive at the node
    std::vector<int16_t> cost(Lx*Ly*cool_down, std::numeric_limits<int16_t>::max());
    std::vector<int16_t> scanned(Lx*Ly, -1);

    std::tuple<int, int> output_node = std::make_tuple(x0, y0);
    
    cost[c3d(x0,y0,0)] = 0;
    int t = 0;
    std::list<std::tuple<int,int>> open_nodes;
    std::list<std::tuple<int,int>> new_nodes;
    int number_of_movements = 0;
    
    open_nodes.push_back(std::make_tuple(x0, y0));
    
    //03-04-24 find possible X, Y of dest (could be more than 1)
    //to know if need to move Up or Down
    //suitable only for case that all outs of same dest are on same Y
    int x_out = x0;
    int y_out = y0;
    for (const auto& output : OutputsWithDest) {
        x_out = std::get<0>(output);
        y_out = std::get<1>(output);
        int z_temp = std::get<2>(output);

        if (z_temp == dest)
            break;
    }
    
    // Forward iterations
    bool found = false;
    while(not found) {
        _total_nodes += (long)open_nodes.size();
        for (const auto& curr_node : open_nodes) {
            int x = std::get<0>(curr_node);
            int y = std::get<1>(curr_node);
            
            //check if node is Output for this dest (suitable also for funnels)
            //if a node is output then destination found and need to stop the while (after finishing this time-step)
            std::tuple<int, int, int> nodeDest = std::make_tuple(x, y, dest);
            auto out = std::find(OutputsWithDest.begin(), OutputsWithDest.end(), nodeDest);
            if (out != OutputsWithDest.end()) {
                found = true;
                output_node = curr_node;
                break; // destination was found
            }
            
            if (t+1 >= cool_down) {
                std::cout << "Panic: path is longer than the value of cool_down = "<< cool_down << " at takt " << t_base << std::endl;
                exit(1);
            }
        
            //03-04-24 split to 2 cases according to the y value of the node compared with the the y value of the out node
            //the difference between the cases is the order of node search
            if (y_out >= y)
            {
                //case 1 - seacrh in the order of up,right,left,stay,down
                if (get_TEG(x, y, t_base + t, dir_up) and cost[c3d(x, y + 1, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x, y + 1, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x, y + 1, t + 1)] = dir_up;
                    if (scanned[c2d(x, y + 1)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y + 1));
                        scanned[c2d(x, y + 1)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_right) and cost[c3d(x + 1, y, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x + 1, y, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x + 1, y, t + 1)] = dir_right;
                    if (scanned[c2d(x + 1, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x + 1, y));
                        scanned[c2d(x + 1, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_left) and cost[c3d(x - 1, y, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x - 1, y, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x - 1, y, t + 1)] = dir_left;
                    if (scanned[c2d(x - 1, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x - 1, y));
                        scanned[c2d(x - 1, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_stay) and cost[c3d(x, y, t + 1)] > cost[c3d(x, y, t)]) {
                    cost[c3d(x, y, t + 1)] = cost[c3d(x, y, t)];
                    from_node[c3d(x, y, t + 1)] = dir_stay;
                    if (scanned[c2d(x, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y));
                        scanned[c2d(x, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_down) and cost[c3d(x, y - 1, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x, y - 1, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x, y - 1, t + 1)] = dir_down;
                    if (scanned[c2d(x, y - 1)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y - 1));
                        scanned[c2d(x, y - 1)] = t;
                    }
                }
            }
            else
            {
                //case 2 - seacrh in the order of down,right,left,stay,up
                if (get_TEG(x, y, t_base + t, dir_down) and cost[c3d(x, y - 1, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x, y - 1, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x, y - 1, t + 1)] = dir_down;
                    if (scanned[c2d(x, y - 1)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y - 1));
                        scanned[c2d(x, y - 1)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_right) and cost[c3d(x + 1, y, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x + 1, y, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x + 1, y, t + 1)] = dir_right;
                    if (scanned[c2d(x + 1, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x + 1, y));
                        scanned[c2d(x + 1, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_left) and cost[c3d(x - 1, y, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x - 1, y, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x - 1, y, t + 1)] = dir_left;
                    if (scanned[c2d(x - 1, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x - 1, y));
                        scanned[c2d(x - 1, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_stay) and cost[c3d(x, y, t + 1)] > cost[c3d(x, y, t)]) {
                    cost[c3d(x, y, t + 1)] = cost[c3d(x, y, t)];
                    from_node[c3d(x, y, t + 1)] = dir_stay;
                    if (scanned[c2d(x, y)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y));
                        scanned[c2d(x, y)] = t;
                    }
                }

                if (get_TEG(x, y, t_base + t, dir_up) and cost[c3d(x, y + 1, t + 1)] > cost[c3d(x, y, t)] + 1) {
                    cost[c3d(x, y + 1, t + 1)] = cost[c3d(x, y, t)] + 1;
                    from_node[c3d(x, y + 1, t + 1)] = dir_up;
                    if (scanned[c2d(x, y + 1)] < t) {
                        new_nodes.push_back(std::make_tuple(x, y + 1));
                        scanned[c2d(x, y + 1)] = t;
                    }
                }
            }
        }
        t++;
        open_nodes = new_nodes;
        new_nodes.clear();
    }

    // Backward iteration and update TEG
    std::list<std::tuple<int,int,int>> path;  // (x,y, dir)
    int x1 = std::get<0>(output_node);
    int y1 = std::get<1>(output_node);
    int x = x1;
    int y = y1;
    int orig_x;
    int orig_y;
    t--;
    int t_out = t_base + t; //save for funnel move
    number_of_movements = int(cost[c3d(x1,y1,t)]);
            
    //loop on all nodes found and construct path for animation + delete from TEG all arcs eliminted by selected path
    while (t>0) {
        orig_x = x;
        orig_y = y;
        switch (from_node[c3d(x,y,t)]) {
            case dir_right:
                orig_x--;
                break;
            case dir_up:
                orig_y--;
                break;
            case dir_left:
                orig_x++;
                break;
            case dir_down:
                orig_y++;
                break;
            default:
                break;
        }
        
        // delete all arcs that leave from this arc's origin node
        for (int d=0; d<5; d++)
            set_TEG(orig_x,orig_y,t_base+t-1,d,0);
        
        // delete all arcs that enter into this arc's destination node
        if (x>0) set_TEG(x-1,y,t_base+t-1,dir_right,0);
        if (x<Lx-1) set_TEG(x+1,y,t_base+t-1,dir_left,0);
        if (y>0) set_TEG(x,y-1,t_base+t-1,dir_up,0);
        if (y<Ly-1) set_TEG(x,y+1,t_base+t-1,dir_down,0);
        set_TEG(x,y,t_base+t-1,dir_stay,0);
        
        // Delete all arcs that enter into this arc's origin node from a direction != current movement direction
        char cur_move_dir = from_node[c3d(x,y,t)];
        if (orig_x > 0 and cur_move_dir != dir_right) set_TEG(orig_x-1,orig_y,t_base+t-1,dir_right,0);
        if (orig_x< Lx-1 and cur_move_dir != dir_left) set_TEG(orig_x+1,orig_y,t_base+t-1,dir_left,0);
        if(orig_y > 0 and cur_move_dir != dir_up) set_TEG(orig_x,orig_y-1,t_base+t-1,dir_up, 0);
        if(orig_y < Ly-1 and cur_move_dir != dir_down) set_TEG(orig_x, orig_y+1, t_base+t-1, dir_down, 0);
        
        // Delete all arcs that leave from this arc's destination node at time t-1 at a direction != current movement direction
        if (cur_move_dir != dir_right) set_TEG(x,y,t_base+t-1,dir_right,0);
        if (cur_move_dir != dir_left) set_TEG(x,y,t_base+t-1,dir_left,0);
        if (cur_move_dir != dir_up) set_TEG(x,y,t_base+t-1,dir_up, 0);
        if (cur_move_dir != dir_down) set_TEG(x,y, t_base+t-1, dir_down, 0);
        
        
        //as everything is stored in the animation file we actual do not need to construct path - we leave it here for now only for QA purposes       
        path.push_front(std::make_tuple(x, y, from_node[c3d(x,y,t)]));
        
        //store moves to animation file
        if (do_animation)
            AnimationData[t_base+t-1].emplace_back(std::make_tuple(count_entered, std::make_tuple(orig_x, orig_y),std::make_tuple(x, y), (x==x1 and y == y1)));
        
        t--;
        x = orig_x;
        y = orig_y;
    }
    
    path.push_front(std::make_tuple(x, y, -1));
    
    if (verbose) {
        std::cout << "(" << x0 << "," << y0 << ") => (" << x1 << "," << y1 << ") : ";
        for (const auto& step : path)
            std::cout << "(" << std::get<0>(step) <<"," << std::get<1>(step) <<" -- " << std::get<2>(step) << ") -> ";
        
        std::cout << std::endl;
        
    }
        
    if (x != x0 or y != y0) {
        std::cout << "Panic: invalid path at takt " << t_base << std::endl;
        exit(1);
    }
    flow_time = int(path.size()-1);
    
    return number_of_movements;
}


int main(int argc, const char * argv[]) {
    
    //get parameters ftom command arguments
    program.add_argument("-V","--verbose")
      .help("increase output verbosity")
      .default_value(false)
      .implicit_value(true);
    
    program.add_argument("-a","--animation")
      .help("Export an animation file")
      .default_value(false)
      .implicit_value(true);

    program.add_argument("-x", "--Lx")
        .help("Horizontal deiminsion of the Grid")
        .default_value(10) // default value if argument is not provided
        .scan<'i', int>(); // 'i' indicates an integer type

    program.add_argument("-y", "--Ly")
        .help("Verical deiminsion of the Grid")
        .default_value(10) // default value if argument is not provided
        .scan<'i', int>(); // 'i' indicates an integer type
       
    program.add_argument("-T","--sim_time")
           .help("Active simulation time not including warmup and cool down. Assumed to be an integer multipliction of --block length")
           .default_value(1000) // default value if argument is not provided
           .scan<'i', int>(); // 'i' indicates an integer type
    
    program.add_argument("-c","--cool_down")
           .help("Number of time steps in cool down period. If not define will use 2*(Lx+Ly)")
           .default_value(0) // default will be calculated as 2*(Lx+Ly)
           .scan<'i', int>(); // 'i' indicates an integer type

    program.add_argument("-s","--seed")
           .help("Random seed (must be non negative)")
           .default_value(0) 
           .scan<'i', int>(); // 'i' indicates an integer type
    
    program.add_argument("-w","--warmup")
           .help("Number of warmup takts (if zero - use 2*(Lx+Ly)")
           .default_value(0) // default will be calculated as 2*(Lx+Ly)
           .scan<'i', int>(); // 'i' indicates an integer type
    
    program.add_argument("-r", "--results")
           .help("Specify the results file")
           .default_value(std::string("results.csv"));

     program.add_argument("-d", "--design")
        .help("Specify the desgin file") //design file includes X,Y dimensions + Input cells + Output cells + Inactive cells
        .default_value(std::string("Grid_Layout_11x11_InS_OutENW.txt"));

    program.add_argument("-H","--header")
      .help("Print header line to the out csv file")
      .default_value(false)
      .implicit_value(true);
  
    program.add_argument("-bl", "--block")
        .help("Block Length for statistical inference")
        .default_value(100)
        .scan<'i', int>(); // 'i' indicates an integer type

 
    try {
      program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << program;
      std::exit(1);
    }

    if (program["--verbose"] == true) {
      std::cout << "Verbosity enabled" << std::endl;
        verbose = true;
    }
    
    Lx = program.get<int>("--Lx");
    Ly = program.get<int>("--Ly");
    T = program.get<int>("--sim_time");
    cool_down = program.get<int>("--cool_down");
    int seed = program.get<int>("--seed");
    int warmup = program.get<int>("--warmup");

    int blockLength = program.get<int>("--block");
    int blockNum = T / blockLength;

    std::list<int> itemDest;
           
    std::ofstream outFile(program.get<std::string>("--results"), std::ios::app);
    //std::ofstream outFileWip(program.get<std::string>("--wipResults"), std::ios::app);
    
    if (outFile.is_open()) {
        if (program["--header"] == true)
            outFile  << std::endl << "file name, Lx x Ly, cells, T, seed, warmup, cool down, number entered, # movements, total flow time, entrance rate, mean movements, mean flow time, std flow time, # of blocks, std block entrance rate, simulation time, time ratio, total nodes"<< std::endl;
    } else {
        std::cerr << "Error: Unable to open file '" << program.get<std::string>("--results") << "' for appending";
        exit(1);
    }
    outFile.close();

    //read data from input file
    std::string fileName = program.get<std::string>("--design");
    std::ifstream inputFile(program.get<std::string>("--design"));

    if (fileName == "")
    {
        std::cerr << "Error: Missing File Name";
        exit(1);

    }

    // Check if the file is open
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return 1; // Return an error code
    }

    // Read and display the content line by line
    std::list<std::tuple<int, int>> Inputs;
    std::list<std::tuple<int, int>> Outputs;
    std::list<std::tuple<int, int, int>> OutputsWithDest;
    std::list<int> Dests;
    std::list<std::tuple<int, int>> Inactive;

    int lineNum = 1;
    std::string line;
    while (std::getline(inputFile, line) and lineNum < 5) {
        if (lineNum == 1)
        {
            //First row of design file contains x,y
            std::string delimiter = ",";

            std::vector<std::string> tokens = splitString(line, delimiter);

            // Display the tokens
            Lx= std::stoi(tokens.at(0));
            Ly = std::stoi(tokens.at(1));
        }
        else if (lineNum == 2)
        {
            //Second row of design file contains input cells
            //Input cell is defined by x,y coordinates
            std::string delimiter = "),(";
            std::vector<std::string> tokens = splitString(line, delimiter);
            for (const auto& token : tokens) {
                std::string tmp = token;
                //remove '(' if first and ')' if last
                if (tmp[0] == '(') tmp = tmp.substr(1);
                if (tmp[tmp.length() - 1] == ')') tmp = tmp.substr(0, tmp.length() - 1);
                std::vector<std::string> tokens2 = splitString(tmp, ",");
                Inputs.push_back(std::make_tuple(std::stoi(tokens2.at(0)), std::stoi(tokens2.at(1))));
            }
        }
        else if (lineNum == 3)
        {
            //Third row of design file contains output cells
            //Output cell is defined by x,y coordinates + Item Type definition
            std::string delimiter = "),(";

            std::vector<std::string> tokens = splitString(line, delimiter);
            for (const auto& token : tokens) {
                std::string tmp = token;
                if (tmp[0] == '(') tmp = tmp.substr(1);
                if (tmp[tmp.length() - 1] == ')') tmp = tmp.substr(0, tmp.length() - 1);
                std::vector<std::string> tokens2 = splitString(tmp, ",");
                Outputs.push_back(std::make_tuple(std::stoi(tokens2.at(0)), std::stoi(tokens2.at(1))));
                OutputsWithDest.push_back(std::make_tuple(std::stoi(tokens2.at(0)), std::stoi(tokens2.at(1)), std::stoi(tokens2.at(2))));
                
                //list of distinct destinations
                auto it = std::find(Dests.begin(), Dests.end(), std::stoi(tokens2.at(2)));
                if (it == Dests.end())
                    Dests.push_back(std::stoi(tokens2.at(2)));
            }
        }
        else if (lineNum == 4 and line.length() > 2)
        {
            //Fourth row of design file contains inactive cells
            //Inactive cell is defined by x,y coordinates
            std::string delimiter = "),(";

            std::vector<std::string> tokens = splitString(line, delimiter);
            for (const auto& token : tokens) {
                std::string tmp = token;
                if (tmp[0] == '(') tmp = tmp.substr(1);
                if (tmp[tmp.length() - 1] == ')') tmp = tmp.substr(0, tmp.length() - 1);
                std::vector<std::string> tokens2 = splitString(tmp, ",");
                //std::cout << "---" << tokens2.at(0) << " " <<  tokens2.at(1);
                Inactive.push_back(std::make_tuple(std::stoi(tokens2.at(0)), std::stoi(tokens2.at(1))));
            }
        }
        lineNum++;
        std::cout << line << std::endl;
    }

    // Close the file
    inputFile.close();

    std::mt19937 rng(seed); // Mersenne Twister generator, seeded with our determinstic seed to allow consistencey
    std::uniform_int_distribution<> rand_dest(0, int(Dests.size()-1));

    if (cool_down == 0)  cool_down = 2 * (Lx + Ly);
    if (warmup == 0) warmup = 2 * (Lx + Ly);
    T1 = T + cool_down + warmup;

    if (program["--animation"] == true) {
        do_animation = true;
        AnimationData.resize(T1);
    }

    std::cout << "Initializing TEG...\n";

    auto start = std::chrono::high_resolution_clock::now();
    TEG.resize(Lx*Ly*T1*5);
    //std::fill_n(std::back_inserter(TEG), Lx*Ly*T1*5, 0); // Fill with zeors - slower
    
    //build TEG from design file data
    for(int x = 0; x <Lx; x++) for(int y=0; y <Ly; y++) for(int t = 0; t < T1; t++) 
    {
        std::tuple<int, int> node = std::make_tuple(x, y);
        //check if node is Output/ Inactive
        auto out = std::find(Outputs.begin(), Outputs.end(), node);
        auto inact = std::find(Inactive.begin(), Inactive.end(), node);
        //block all arcs if current node is output or inactive
        if (out != Outputs.end() || inact != Inactive.end())
        {
            set_TEG(x, y, t, dir_up, 0);
            set_TEG(x, y, t, dir_right, 0);
            set_TEG(x, y, t, dir_left, 0);
            set_TEG(x, y, t, dir_down, 0);
            set_TEG(x, y, t, dir_stay, 0);
        }
        else
        {
            //add arcs for stay
            set_TEG(x, y, t, dir_stay, 1);
            //add arcs for up move if next_node not in Input or Inactive
            if (y < Ly - 1)
            {
                std::tuple<int, int> next_node = std::make_tuple(x, y + 1);
                auto in = std::find(Inputs.begin(), Inputs.end(), next_node);
                auto inact2 = std::find(Inactive.begin(), Inactive.end(), next_node);
                if (in == Inputs.end() and inact2 == Inactive.end())
                    set_TEG(x, y, t, dir_up, 1);
                else set_TEG(x, y, t, dir_up, 0);
            }
            else set_TEG(x, y, t, dir_up, 0);

            //add arcs for right move if next_node not in Input or Inactive
            if (x < Lx - 1)
            {
                std::tuple<int, int> next_node = std::make_tuple(x + 1, y);
                auto in = std::find(Inputs.begin(), Inputs.end(), next_node);
                auto inact2 = std::find(Inactive.begin(), Inactive.end(), next_node);
                if (in == Inputs.end() and inact2 == Inactive.end())
                    set_TEG(x, y, t, dir_right, 1);
                else set_TEG(x, y, t, dir_right, 0);
            }
            else set_TEG(x, y, t, dir_right, 0);

            //add arcs for left move if next_node not in Input or Inactive
            if (x > 0)
            {
                std::tuple<int, int> next_node = std::make_tuple(x - 1, y);
                auto in = std::find(Inputs.begin(), Inputs.end(), next_node);
                auto inact2 = std::find(Inactive.begin(), Inactive.end(), next_node);
                if (in == Inputs.end() and inact2 == Inactive.end())
                    set_TEG(x, y, t, dir_left, 1);
                else set_TEG(x, y, t, dir_left, 0);
            }
            else set_TEG(x, y, t, dir_left, 0);

            //add arcs for down move if next_node not in Input or Inactive
            if (y > 0)
            {
                std::tuple<int, int> next_node = std::make_tuple(x, y - 1);
                auto in = std::find(Inputs.begin(), Inputs.end(), next_node);
                auto inact2 = std::find(Inactive.begin(), Inactive.end(), next_node);
                if (in == Inputs.end() and inact2 == Inactive.end())
                    set_TEG(x, y, t, dir_down, 1);
                else set_TEG(x, y, t, dir_down, 0);
            }
            else set_TEG(x, y, t, dir_down, 0);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Time taken to initialize TEG: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "Simulating rectangular " << Lx << "x" << Ly << " system with " << Inputs.size() << " inputs and " << Outputs.size() << " outputs" << "  T=" << T << "  warmup=" << warmup << "  cool_down=" << cool_down <<std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    int total_number_of_movements = 0;
    int total_flow_time = 0;
    double total_flow_time_square = 0;
    
    int flow_time = 0;
    int count_before_warmup = 0;
    long int total_nodes = 0; // total number of scanned nodes
    std::vector<double_t> items_enter(blockNum, 0);
    
    //Simulating dyanmic sorting
    //loop on t to enter new units into the grid
    for (int t = 0; t < T+warmup; t++) {
        if (t % ((T + warmup)/10) == 0 and t > 0) std::cout << 100 * t / (T + warmup) << "% done." << std::endl;
        if (verbose)  std::cout << "************ Takt:" << t << " ************** " << std::endl;
        if (t==warmup) {
            count_before_warmup = count_entered;
            total_number_of_movements = 0;
            total_flow_time = 0;
            total_flow_time_square = 0;
        }
        //loop on Inputs to find an open input cell to enter a new unit
        for (const auto& input : Inputs) {
            int x0 = std::get<0>(input);
            int y0 = std::get<1>(input);
            if (get_TEG(x0,y0,t, dir_stay)) {
                
                //for each open input cell randomaly selection a destination and find shortest path on TEG
                int dest = rand_dest(rng);
                
                if (do_animation) itemDest.emplace_back(dest);
                total_number_of_movements += find_path(t, x0, y0 , dest, OutputsWithDest, flow_time, total_nodes);
                total_flow_time += flow_time;

                if (t >= warmup) items_enter[(t-warmup)/blockLength] += 1; //save items entered per block
                total_flow_time_square += std::pow(flow_time, 2); 
                count_entered++;
            }
        }
    }
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    int entered_after_warmup = count_entered- count_before_warmup;
    std::cout << "Simulation time: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "Total number of items entered: " << count_entered << " in " << T+warmup << " takts. The rate is " << float(count_entered) / float(T + warmup) << std::endl;
   
    std::cout << "Total number of items entered after warmup: " << entered_after_warmup << " in " << T << " takts. The rate (after warmup) is " << float(entered_after_warmup) / float(T) << std::endl;
    
    std::cout << "Sim time / realtime ratio (takt = second): " << float(duration.count()) / (float((T + warmup)*1000)) << " in " << T + warmup << " takts (less than 1 is good)" << std::endl;
    
    std::cout << "Total number of movements (of items that entered aftetr warmup): " << total_number_of_movements << std::endl;
    std::cout << "Total flow time (of items enter aftetr warmup): " << total_flow_time << std::endl;
    std::cout << "Mean flow time: " << float(total_flow_time) / entered_after_warmup << std::endl;
    std::cout << "Total number of nodes: " << total_nodes << std::endl;
    
    double mean = 0;
    int L = count_entered - count_before_warmup;
    double stdBlockEntry = calculateStandardDeviation(items_enter, mean);
    double stdFlowTime = std::sqrt(total_flow_time_square / L - std::pow(flow_time / L, 2));
   
    //save data to results file
    //"file name, Lx x Ly, cells, T, seed, warmup, cool down, number entered, # movements, total flow time, entrance rate, mean movements, mean flow time, std flow time, # of blocks , std block entrance rate, simulation time, time ratio, total nodes"
    outFile.open(program.get<std::string>("--results"), std::ios::app);
    outFile << fileName << "," << Lx <<"x" << Ly << "," << Lx*Ly << "," << T << "," << seed << "," << warmup << "," << cool_down << "," << L << "," << total_number_of_movements << "," << total_flow_time << "," << float(entered_after_warmup) / float(T) << "," << float(total_number_of_movements) / entered_after_warmup << "," << float(total_flow_time) / entered_after_warmup << "," << stdFlowTime << "," << items_enter.size() << "," << stdBlockEntry << "," << duration.count() << "," << float(duration.count()) / (T*1000) << "," << total_nodes << std::endl;

    outFile.close();
    
    //save all paths to script file for animation
    if (do_animation) {
        std::string script_file_name("script_"+std::to_string(Lx)+"x"+std::to_string(Ly)+"_T"+std::to_string(T) + ".txt");
        
        std::ofstream scriptFile(script_file_name);
        
        if (scriptFile.is_open()) {
            scriptFile << Lx << std::endl << Ly << std::endl;
            
            // Export I
            int index = 0;
            scriptFile << "[";
            for (const auto& input : Inputs) {
                int x = std::get<0>(input);
                int y = std::get<1>(input);

                if (index > 0) 
                    scriptFile << ",";
                scriptFile << "(" << x << "," << y << ")";
                index++;
            }
            scriptFile << "]" << std::endl;
            
            // export OT
            index = 0;
            scriptFile << "[";
            for (const auto& output : OutputsWithDest) {
                int x = std::get<0>(output);
                int y = std::get<1>(output);
                int z = std::get<2>(output);

                if (index > 0)
                    scriptFile << ",";
                scriptFile << "(" << x << ","<< y << "," << z <<")";
                index++;
            }
            scriptFile << "]" << std::endl;
            
            // export IA. Blocked locations
            index = 0;
            scriptFile << "[";
            for (const auto& inactive : Inactive) {
                int x = std::get<0>(inactive);
                int y = std::get<1>(inactive);

                if (index > 0)
                    scriptFile << ",";
                scriptFile << "(" << x << "," << y << ")";
                index++;
            }
            scriptFile << "]" << std::endl;

            scriptFile << count_entered << std::endl; // export count
            
            scriptFile << "[";   // export script
            for (int t = 0; t < T1; ++t) {
                if (AnimationData[t].size() >0) {
                    if (t!=0) scriptFile << ",";
                    scriptFile << "[";
                    bool lFirst = true;
                    for( auto &takt_data : AnimationData[t]){
                        if (lFirst) lFirst = false;
                        else scriptFile << ",";
                        std::string lastMove = (std::get<3>(takt_data)) ? "True" :"False";
                        scriptFile << "(" << std::get<0>(takt_data) << ",(" << std::get<0>(std::get<1>(takt_data)) << "," << std::get<1>(std::get<1>(takt_data)) << "),(" << std::get<0>(std::get<2>(takt_data)) << "," << std::get<1>(std::get<2>(takt_data)) << ")," << lastMove << ")";
                    }
                    scriptFile << "]";
                }
            }
            scriptFile << "]" << std::endl;
            
            scriptFile << "{";
            int count = 0;
            for( auto &item : itemDest) {
                if (count!=0) scriptFile << ", ";
                scriptFile << count << ": " << item;
                count++;
            }
            scriptFile << "}" << std::endl;
        } else {
            std::cout << "Panic: can't open script file '" << script_file_name << "' for writing" << std::endl;
            exit(1);
        }
    }
    TEG.clear();
    return 0;
}
