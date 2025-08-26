# help functions for mapping metadata to masterdata object types
def metadata_to_instance(metadata: dict, instance):
    props = metadata_to_masterdata(metadata, instance)
    for k, v in props.items():
        setattr(instance, k, v)
    return instance


def metadata_to_masterdata(metadata: dict, object_instance):
    avaliable_values = object_instance._property_metadata
    data = metadata
    object_prop_list = {
        prop_name: data_value
        for prop_name, prop_obj in avaliable_values.items()
        if (data_value := data.get(prop_obj.code)) is not None
    }
    return object_prop_list
