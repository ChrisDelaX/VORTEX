function[cellH]=CreateCell(xCell,yCell,cellColor,nbQuadrants)

lwz = 2;
line(xCell,yCell,'color',cellColor,'linewidth',lwz)
if nbQuadrants >= 2, line(-xCell,yCell,'color',cellColor,'linewidth',lwz)
end
if nbQuadrants == 4,
    line(-xCell,-yCell,'color',cellColor,'linewidth',lwz)
    line(xCell,-yCell,'color',cellColor,'linewidth',lwz)
end




