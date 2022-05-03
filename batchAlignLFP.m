% batch process align LFP

% running analysis for each alingment. 
for as = 1:2
    switch as
        case 1
            % choosing alignment spot: stimulus or resoponse
            alignSpot = 'stimulus';
        case 2
            % choosing alignment spot: stimulus or resoponse
            alignSpot = 'response';
    end
    % looping over patients:
    for pt = 1:2
        switch pt
            case 1
                % YDX1
                ptID = 'YDX1';
                parentDir = 'D:\Data\Elliot\AlgoPlaceCells\NotBirds_data';
                run1_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-004_task-Birds-Object-run-01_blk-01\EMU-004_subj-YDX_task-Birds-Object_run-01_blk-01_NSP-1.ns3';
                run2_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-005_task-Birds-Object-run-01_blk-02\EMU-005_subj-YDX_task-Birds-Object_run-01_blk-02_NSP-1.ns3';

            case 2
                % YEA1
                ptID = 'YEA1';
                parentDir = 'D:\Data\Elliot\AlgoPlaceCells\NotBirds_data';
                run1_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-011_task-Birds-Object-run-01_blk-01\EMU-011_task-Birds-Object-run-01_blk-01\EMU-011_subj-YEA_task-Birds-Object-run-01_blk-01_NSP-1.ns3';
                run2_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-012_task-Birds-Object-run-01_blk-02\EMU-012_task-Birds-Object-run-01_blk-02\EMU-012_subj-YEA_task-Birds-Object-run-01_blk-02_NSP-1.ns3';
        end

        % aligning LFP
        alignLFPandSpikes(parentDir,ptID,alignSpot,run1_fName,run2_fName);

    end
end