from crams.crams import (
    FragmentationModel_Fluka4Dragon,
    InelasticModel_Glauber,
    Input,
    Runner,
    git_sha1,
    ParticleList,
)

print(git_sha1())

runner = Runner(InelasticModel_Glauber, FragmentationModel_Fluka4Dragon, ParticleList())

params = Input(H_kpc=3)
print(params.describe())
result = runner.compute(params, dumpToFile=False, verbose=True, ignoreInputInitParams=True)

print(result)
