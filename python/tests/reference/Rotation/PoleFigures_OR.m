% Start MTEX first in Matlab

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

symmetry = {crystalSymmetry('m-3m', [1 1 1], 'mineral', 'Iron', 'color', 'light blue')}

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');


lattice_types = {'BCC','FCC'};
models        = {'Bain','GT','GT_prime','KS','NW','Pitsch'};

rotation = containers.Map;
rotation('BCC') = 'Passive Rotation';
rotation('FCC') = 'Active Rotation';

for lattice = lattice_types
  for p = 0:length(models)/3-1
    EBSD_data = {loadEBSD(strcat(lattice,'_',models{p*3+1},'.txt'),symmetry,'interface','generic',...
                          'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', rotation(char(lattice))),
                 loadEBSD(strcat(lattice,'_',models{p*3+2},'.txt'),symmetry,'interface','generic',...
                          'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', rotation(char(lattice))),
                 loadEBSD(strcat(lattice,'_',models{p*3+3},'.txt'),symmetry,'interface','generic',...
                          'ColumnNames', { 'phi1' 'Phi' 'phi2' 'x' 'y'}, 'Bunge', rotation(char(lattice)))}
    h = [Miller(1,0,0,symmetry{1}),Miller(1,1,0,symmetry{1}),Miller(1,1,1,symmetry{1})];            % 3 pole figures
    plotPDF(EBSD_data{1}.orientations,h,'MarkerSize',5,'MarkerColor','r','DisplayName',models{p*3+1})
    hold on
    plotPDF(EBSD_data{2}.orientations,h,'MarkerSize',4,'MarkerColor','b','DisplayName',models{p*3+2})
    plotPDF(EBSD_data{3}.orientations,h,'MarkerSize',3,'MarkerColor','g','DisplayName',models{p*3+3})
    legend('show','location','southoutside','Interpreter', 'none') 
    orient('landscape')
    print('-bestfit',strcat(int2str(p+1),'_',char(lattice),'.pdf'),'-dpdf')
    close
  end
end