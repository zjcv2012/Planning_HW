function [startQ goalQ] = getValidConf(mapfile);
    LINKLENGTH_CELLS=10;
    envmap = load(mapfile);

    close all;
    
    while(1)
        [startQ goalQ] = generateRandomConf();
        armplan = armplanner(envmap, startQ, goalQ, 0);
        if(size(armplan,1)>1)
            break;
        end
    end
end