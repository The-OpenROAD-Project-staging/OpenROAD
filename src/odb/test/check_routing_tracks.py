import openroad as ord
import helpers
import odb

tech = ord.Tech()
tech.readLEF("data/Nangate45/NangateOpenCellLibrary.mod.lef")
db = ord.get_db()
odb.read_def(db, "data/gcd/floorplan.def")

# This has layer info etc
dbtech = odb.dbTech()

chip = db.getChip()
block = chip.getBlock()


for i in range(0, dbtech.getRoutingLayerCount() ):
    routingLayer = dbtech.findRountingLayer(i)

    if routingLayer == None:
        print("FAIL: Read routing layer Failed")
        exit(1)

    routingTrack = block.findTrackGrid(routingLayer)

    if routingTrack == None:
        print("FAIL: Read routing track for layer $l Failed")
        exit(1)
        
print("pass")
exit(0)
