#include "run3ScoutingLooper.C"

unsigned int numberOfFilesToRun = 1;


int main(int argc, char **argv) {
  // Sample list: Data
  std::vector<TString> dataFiles;
  const fs::path dataDir{"/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/"};
  unsigned int iFile=0;
  for (auto const& dataFile : fs::directory_iterator{dataDir}) {
    if (!TString(dataFile.path()).Contains(".root"))
      continue;
    else
      dataFiles.push_back(TString(dataFile.path()));
    iFile++;
    if (iFile == numberOfFilesToRun)
      break;
  }

  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << dataFiles.size() << "\n";
  std::cout << "################################## \n";
  run3ScoutingLooper(dataFiles, "2022", "Data");
  //run3ScoutingLooper({"/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/output_340.root"}, "2022", "Data", "temp_data");
  //run3ScoutingLooper({"/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/output_10.root","/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/output_102.root"}, "2022", "Data", "temp_data");

  return 0;
}

