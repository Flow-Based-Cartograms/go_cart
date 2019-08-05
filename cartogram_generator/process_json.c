/****************************** Inclusions. **********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <cjson/cJSON.h>
#include "cartogram.h"

/**************************** Function prototypes. ***************************/

void process_json(char *json_file_name, char *program_name);

static int compare(const void* a, const void* b) 
  { 
  const name_id_pair *da = a, *db = b;
  return strcmp(da->name, db->name);
} 

void process_json(char *json_file_name, char *program_name){
	if(json_file_name == NULL){
      fprintf(stderr, "ERROR: Cannot find GeoJSON file.\n");
      exit(1);
    }else{
      FILE *json_file;
      char tmp_file_name[100];
      char extracted_file_name[100];
      if ((json_file = fopen(json_file_name,"r")) == NULL) {
        fprintf(stderr,"ERROR: Cannot find GeoJSON file.\n");
        exit(1);
      }else{
        size_t namelen = strlen(json_file_name), extlenjson = strlen(".json"), extlengeojson = strlen(".geojson");
        if(namelen >= extlenjson && !strcmp(json_file_name + namelen - extlenjson, ".json")){
          strncpy(tmp_file_name, json_file_name, sizeof(tmp_file_name) - strlen(tmp_file_name) - 1);
          strncpy(extracted_file_name, tmp_file_name, namelen - extlenjson);
          extracted_file_name[namelen - extlenjson] = '\0';
        }else if(namelen >= extlengeojson && !strcmp(json_file_name + namelen - extlengeojson, ".geojson")){
          strncpy(tmp_file_name, json_file_name, sizeof(tmp_file_name) - strlen(tmp_file_name) - 1);
          strncpy(extracted_file_name, tmp_file_name, namelen - extlengeojson);
          extracted_file_name[namelen - extlengeojson] = '\0';
        }else{
          fprintf(stderr,"ERROR: Map file not in proper file format. Map file needs to be a .json file or .geojson file.\n");
          exit(1);
        }
        fprintf(stderr, "Processing GeoJSON File...\n");
        char *strJson = NULL;
        int line_len;
        int MAX_STRING_LENGTH_JSON = 1000000;
        char* str = (char *) malloc(MAX_STRING_LENGTH_JSON);
        int strJson_size = 0;

        while (fgets(str, MAX_STRING_LENGTH_JSON, json_file) != NULL){
          line_len = strlen(str);
          if(line_len > 0 && str[line_len - 1] == '\n'){
            str[line_len - 1] = '\0';
          }

          if(strJson == NULL){
            strJson_size = line_len + 1;
            strJson = (char *) malloc(line_len + 1);
            strncpy(strJson, str, strJson_size - strlen(strJson) - 1);
          }else{
            strJson_size = (strlen(strJson)) + line_len + 1;
            strJson = (char *) realloc(strJson, (strlen(strJson)) + line_len + 1);
            strncat(strJson, str, strJson_size - strlen(strJson) - 1);
          }
        }
        free(str);
        cJSON *root = cJSON_Parse(strJson);
        free(strJson);
        cJSON *features = cJSON_GetObjectItemCaseSensitive(root, "features");
        cJSON *bbox = cJSON_GetObjectItemCaseSensitive(root, "bbox");
        n_reg = cJSON_GetArraySize(features);
        
        cJSON *feature0 = cJSON_GetArrayItem(features, 0);
        cJSON *feature1 = cJSON_GetArrayItem(features, 1);
        cJSON *feature2 = cJSON_GetArrayItem(features, 2);
        cJSON *feature_properties_0 = cJSON_GetObjectItemCaseSensitive(feature0, "properties");
        cJSON *feature_properties_1 = cJSON_GetObjectItemCaseSensitive(feature1, "properties");
        cJSON *feature_properties_2 = cJSON_GetObjectItemCaseSensitive(feature2, "properties");
        char property_name[50] = "NAME_1";
        cJSON *feature_name_obj = cJSON_GetObjectItemCaseSensitive(feature_properties_0, property_name);
        cJSON *feature_country_gid = cJSON_GetObjectItemCaseSensitive(feature_properties_0, "GID_0");
        if(!feature_name_obj){
          fprintf(stderr, "\n[Property]: [Values]\n\n");
          cJSON *property_iterator = NULL;
          cJSON_ArrayForEach(property_iterator, feature_properties_0){
            char *key = property_iterator->string;
            char *value0 = property_iterator->valuestring;
            char *value1 = cJSON_GetObjectItemCaseSensitive(feature_properties_1, key)->valuestring;
            char *value2 = cJSON_GetObjectItemCaseSensitive(feature_properties_2, key)->valuestring;
            fprintf(stderr, "%s: %s, %s, %s, ...\n", key, value0, value1, value2);
          }
          fprintf(stderr,"\nPlease identify the property which contains the names of each region: ");
          if(scanf("%s", property_name) != 1){
            fprintf(stderr, "Failed to read input.\n");
            exit(1);
          }
          feature_name_obj = cJSON_GetObjectItemCaseSensitive(feature_properties_0, property_name);
          while(!feature_name_obj){
            fprintf(stderr,"Input does not match to any property. Please try again: ");
            if(scanf("%s", property_name) != 1){
              fprintf(stderr, "Failed to read input.\n");
              exit(1);
            }
            feature_name_obj = cJSON_GetObjectItemCaseSensitive(feature_properties_0, property_name);
          }
        }
        cJSON *feature_iterator = NULL;
        name_id_pair feature_name_pair[n_reg];
        int ctr = 0;
        cJSON_ArrayForEach(feature_iterator, features){
          cJSON *feature_properties = cJSON_GetObjectItemCaseSensitive(feature_iterator, "properties");
          cJSON *feature_nameproperty = cJSON_GetObjectItemCaseSensitive(feature_properties, property_name);
          strcpy(feature_name_pair[ctr].name , feature_nameproperty->valuestring);
          feature_name_pair[ctr].id = ctr;
          ctr++;
        }
        qsort(feature_name_pair, n_reg, sizeof *feature_name_pair, compare); 
        char *csv_file_name = malloc(sizeof *csv_file_name * 100);
        if(feature_country_gid){
          char *gid = malloc(sizeof *gid * strlen(feature_country_gid->valuestring));
          strcpy(gid, feature_country_gid->valuestring);
          for(int i=0; gid[i]; i++){
            gid[i] = tolower(gid[i]);
          }
          strcpy(extracted_file_name, gid);
        }
        strcpy(csv_file_name, extracted_file_name);
        strcat(csv_file_name, "_data.csv");
        ctr = 1;
        FILE *tmp_file;
        while((tmp_file = fopen(csv_file_name,"r")) != NULL){
          int length = snprintf(NULL, 0, "%d", ctr);
          char* suffix = malloc(length + 1);
          snprintf(suffix, length + 1, "%d", ctr);
          strcpy(csv_file_name, extracted_file_name);
          strcat(csv_file_name, "_data_(");
          strcat(csv_file_name, suffix);
          strcat(csv_file_name, ").csv");
          ctr++;;
        }
        FILE *csv_file = fopen(csv_file_name, "w");
        fprintf(csv_file, "Region Id,Region Data,Region Name\n");
        for (int i=0; i<n_reg; i++){
          cJSON *current_feature = cJSON_GetArrayItem(features, feature_name_pair[i].id);
          cJSON *current_feature_properties = cJSON_GetObjectItemCaseSensitive(current_feature, "properties");
          int length = snprintf(NULL, 0, "%d", i+1);
          char* region_id_str = malloc(length + 1);
          snprintf(region_id_str, length + 1, "%d", i+1);
          if(cJSON_GetObjectItemCaseSensitive(current_feature_properties, "cartogram_id") != NULL){
            cJSON_DetachItemFromObject(current_feature_properties, "cartogram_id");
          }
          cJSON_AddStringToObject(current_feature_properties, "cartogram_id", region_id_str);
          fprintf(csv_file, "%s,,%s\n", region_id_str, feature_name_pair[i].name);
          free(region_id_str);
        }
        fflush(csv_file);
        fclose(csv_file);
        if(!bbox){
          double json_minx = INFINITY;
          double json_maxx = -INFINITY;
          double json_miny = INFINITY;
          double json_maxy = -INFINITY;
          double x;
          double y;
          feature_iterator = NULL;
          cJSON_ArrayForEach(feature_iterator, features){
            cJSON *feature_geometry = cJSON_GetObjectItemCaseSensitive(feature_iterator, "geometry");
            char * type = cJSON_GetObjectItemCaseSensitive(feature_geometry, "type")->valuestring;
            cJSON * feature_geometry_coordinates = cJSON_GetObjectItemCaseSensitive(feature_geometry, "coordinates");
            if(strcmp(type, "Polygon") == 0){
              cJSON * polygon_array_of_linear_rings = feature_geometry_coordinates;
              cJSON * linear_ring_array_of_positions = NULL;
              cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
                cJSON * position = NULL;
                cJSON_ArrayForEach(position, linear_ring_array_of_positions){
                  x = cJSON_GetArrayItem(position, 0)->valuedouble;
                  y = cJSON_GetArrayItem(position, 1)->valuedouble;
                  json_minx = MIN(json_minx, x);
                  json_maxx = MAX(json_maxx, x);
                  json_miny = MIN(json_miny, y);
                  json_maxy = MAX(json_maxy, y);
                }
              }
            }else if(strcmp(type, "MultiPolygon") == 0){
              cJSON * multipolygon_array_of_polygons = feature_geometry_coordinates;
              cJSON * polygon_array_of_linear_rings = NULL;
              cJSON_ArrayForEach(polygon_array_of_linear_rings, multipolygon_array_of_polygons){
                cJSON * linear_ring_array_of_positions = NULL;
                cJSON_ArrayForEach(linear_ring_array_of_positions, polygon_array_of_linear_rings){
                  cJSON * position = NULL;
                  cJSON_ArrayForEach(position, linear_ring_array_of_positions){
                    x = cJSON_GetArrayItem(position, 0)->valuedouble;
                    y = cJSON_GetArrayItem(position, 1)->valuedouble;
                    json_minx = MIN(json_minx, x);
                    json_maxx = MAX(json_maxx, x);
                    json_miny = MIN(json_miny, y);
                    json_maxy = MAX(json_maxy, y);
                  }
                }
              }
            }else{
              fprintf(stderr, "Error: Region contains geometry other than polygons and multipolygons.\n");
              exit(1);
            }
          }
          bbox = cJSON_AddArrayToObject(root, "bbox");
          cJSON_AddItemToArray(bbox, cJSON_CreateNumber(json_minx));
          cJSON_AddItemToArray(bbox, cJSON_CreateNumber(json_miny));
          cJSON_AddItemToArray(bbox, cJSON_CreateNumber(json_maxx));
          cJSON_AddItemToArray(bbox, cJSON_CreateNumber(json_maxy));
        }
        char *processed_json_file_name = malloc(sizeof *processed_json_file_name * 100);
        strcpy(processed_json_file_name, extracted_file_name);
        strcat(processed_json_file_name, "_processedmap.json");
        ctr = 1;
        while((tmp_file = fopen(processed_json_file_name,"r")) != NULL){
          int length = snprintf(NULL, 0, "%d", ctr);
          char* suffix = malloc(length + 1);
          snprintf(suffix, length + 1, "%d", ctr);
          strcpy(processed_json_file_name, extracted_file_name);
          strcat(processed_json_file_name, "_processedmap_(");
          strcat(processed_json_file_name, suffix);
          strcat(processed_json_file_name, ").json");
          ctr++;;
        }
        FILE *processed_json_file = fopen(processed_json_file_name, "w");
        char *json_output = cJSON_Print(root);
        cJSON_Delete(root);
        fprintf(processed_json_file, "%s", json_output);
        fflush(processed_json_file);
        fclose(processed_json_file);
        free(json_output);
        fprintf(stderr, "\nExported map file: %s.\n", processed_json_file_name);
        fprintf(stderr, "Exported area file: %s. Please fill in region data in this file.\n", csv_file_name);
        fprintf(stderr, "Please use both these files for the cartogram generation.\n\n");
        print_usage(program_name);
      }
      return;
    }
}