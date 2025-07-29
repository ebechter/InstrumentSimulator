%SiliconGrating
clear
close

directory = 'S:\Simulator\RefFiles\Curves\RawCurves\echelle80k7_Littrow';
type = '.dat'; % need to change this to a check
fileID = strcat(directory,'/*',type);
allFiles = dir(fileID);
allNames = {allFiles.name};
noFrames = size(allNames,2); %M is the number of files in the directory

figure
hold on
for ii = 1:length(allNames)
    
SiliconR681(:,:,ii) = importGratingFile(allNames{ii});

end

figure()
plot(SiliconR681(:,1,3)*1e9,SiliconR681(:,2,15))
hold on
plot(SiliconR681(:,1,3)*1e9,SiliconR681(:,3,15))

%Two of the orders are not complete! This is why we aren't using this data
%right now.