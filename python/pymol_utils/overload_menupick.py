'''ouptut the atom index when clicking on an atom, instead of the /chainName/resName/atomName
'''
python
import pymol.menu
original_pick_menu = pymol.menu.pick_menu
def custom_pick_menu(cmd, sele1, sele2):
    sele_ids = []
    cmd.iterate(sele2, 'sele_ids[:] = [ID, index]', space=locals())
    sele1 += ' ID:%d index:%d' % tuple(sele_ids)
    return original_pick_menu(cmd, sele1, sele2)
pymol.menu.pick_menu = custom_pick_menu
python end
