from fuzzySys import DiffSystem, MemberFunct
from numpy import inf
#
Ruido_Alto = MemberFunct('Ruido', 'Alto', False, None, 80, 90, inf, inf)
Ruido_Medio = MemberFunct('Ruido', 'Medio', False, None, 50, 60, 80, 90)
Ruido_Bajo = MemberFunct('Ruido', 'Bajo', False, None, 0, 0, 50, 60)

VolumenDeLaSignal_Alto = MemberFunct('VolumenDeLaSignal', 'Alto', False, None, 70, 90, inf, inf)
VolumenDeLaSignal_Medio = MemberFunct('VolumenDeLaSignal', 'Medio', True, None, 50, 70, 90)
VolumenDeLaSignal_Bajo = MemberFunct('VolumenDeLaSignal', 'Bajo', False, None, 0, 0, 50, 60)

VolumenDeLaTV_A = MemberFunct('VolumenDeLaTV','Bajo',False,None,0,0,5,10)
VolumenDeLaTV_B = MemberFunct('VolumenDeLaTV','Medio',True,None,5,10,15)

mapField_MemberToFunc = {
    ('Ruido', 'Alto'): Ruido_Alto,
    ('Ruido', 'Medio'): Ruido_Medio,
    ('Ruido', 'Bajo'): Ruido_Bajo,
    ('VolumenDeLaSignal', 'Alto'): VolumenDeLaSignal_Alto,
    ('VolumenDeLaSignal', 'Medio'): VolumenDeLaSignal_Medio,
    ('VolumenDeLaSignal', 'Bajo'): VolumenDeLaSignal_Bajo,
    ('VolumenDeLaTV','Bajo'): VolumenDeLaTV_A,
    ('VolumenDeLaTV','Medio'): VolumenDeLaTV_B,
}
input_data = {
    'Ruido': 87,
    'VolumenDeLaSignal': 80
}

diffSys = DiffSystem(mapField_MemberToFunc)
diffSys.add_rule(diffSys.fuzzy_or('Ruido Alto','VolumenDeLaSignal Alto'), 'VolumenDeLaTV Bajo')
diffSys.add_rule(diffSys.fuzzy_and('Ruido Alto', diffSys.fuzzy_not('VolumenDeLaSignal Alto')), 'VolumenDeLaTV Medio')
diffSys.add_rule(diffSys.fuzzy_and('Ruido Medio', diffSys.fuzzy_not('VolumenDeLaSignal Alto')), 'VolumenDeLaTV Medio')
diffSys.add_rule(diffSys.fuzzy_and('Ruido Bajo', 'VolumenDeLaSignal Medio'), 'VolumenDeLaTV Bajo')
diffSys.add_rule(diffSys.fuzzy_and('Ruido Bajo', 'VolumenDeLaSignal Bajo'), 'VolumenDeLaTV Medio')

result, rules_evaluations = diffSys.get_output(input_data, 'mamdani', 'mom')
print(result)
print(rules_evaluations)
